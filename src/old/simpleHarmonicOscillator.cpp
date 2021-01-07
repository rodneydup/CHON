#include "al/app/al_App.hpp"
#include "Gamma/Oscillator.h"
#include "al/graphics/al_Shapes.hpp"
#include "al/ui/al_ControlGUI.hpp"
#include "al/ui/al_Parameter.hpp"
#include "al/math/al_Ray.hpp"
#include "al/ui/al_BoundingBox.hpp"
#include "al/ui/al_Pickable.hpp"
#include <vector>


using namespace al;

struct MyApp: public App {
    float color = 0.0;
    Mesh particle;
    Mesh boundaries[2];
    
    Parameter k{"k", "", 1.0f, "", 0.1f, 10.0f}; // spring constant
    Parameter m{"m", "", 10.0f, "", 0.1f, 10.0f}; // mass
    Parameter b{"b", "", 0.0f, "", 0.0f, 10.0f}; // damping

    ControlGUI gui;

    float x = 0.5; // position
    float v; // velocity
    float w; // angular velocity

    float boundaryL = -1;
    float boundaryR = 1;

    float step; // time step

    void onInit() { // Called on app start
        std::cout << "onInit()" << std::endl;
    }

    void onCreate() { // Called when graphics context is available
        std::cout << "onCreate()" << std::endl;
        gui << k << m << b;
        gui.init();
        addIcosphere(particle, 0.1, 4);

        // Create face normals
        particle.decompress();
        particle.generateNormals();

        for (int i = 0; i < 2; i++) {
            addCube(boundaries[i], 0, 0.1);
            boundaries[i].translate(i*2-1, 0, -5);
            boundaries[i].decompress();
            boundaries[i].generateNormals();
        }

        navControl().useMouse(false);
    }

    void onAnimate(double dt) { // Called once before drawing
        float a = ((-k * x) / m) - (v * b); // acceleration
        v += a/60; // update velocity dividing by frame rate
        x += v; // update position
        if (x < boundaryL) x = boundaryL;
        if (x > boundaryR) x = boundaryR;
    } 

    void onDraw(Graphics &g) { // Draw function
        g.clear(0, 0, 0);

        g.depthTesting(true);
        g.lighting(true);

        g.color(0.5, 0.5, 0.5);
        g.polygonMode(Graphics::FILL);
        g.draw(boundaries[0]);
        g.draw(boundaries[1]);
        g.translate(x, 0, -5);

        g.draw(particle);
        gui.draw(g);
    }

    void onSound(AudioIOData &io) { // Audio callback
        while (io()) {
            
        }
    }

    void onMessage(osc::Message &m) { // OSC message callback
        m.print();
    }

    // helper functions to get scene ray from mouse events
  Vec3d unproject(Vec3d screenPos) {
    auto &g = graphics();
    auto mvp = g.projMatrix() * g.viewMatrix() * g.modelMatrix();
    Matrix4d invprojview = Matrix4d::inverse(mvp);
    Vec4d worldPos4 = invprojview.transform(screenPos);
    return worldPos4.sub<3>(0) / worldPos4.w;
  }

  Rayd getPickRay(int screenX, int screenY) {
    Rayd r;
    Vec3d screenPos;
    screenPos.x = (screenX * 1. / width()) * 2. - 1.;
    screenPos.y = ((height() - screenY) * 1. / height()) * 2. - 1.;
    screenPos.z = -1.;
    Vec3d worldPos = unproject(screenPos);
    r.origin().set(worldPos);

    screenPos.z = 1.;
    worldPos = unproject(screenPos);
    r.direction().set(worldPos);
    r.direction() -= r.origin();
    r.direction().normalize();
    return r;
  }
  bool onMouseMove(const Mouse &m) override {
    // make a ray from mouse location
    Rayd r = getPickRay(m.x(), m.y());
    return true;
  }
  bool onMouseDown(const Mouse &m) override {
    Rayd r = getPickRay(m.x(), m.y());
    if (float(m.y())/height() > 0.4 && float(m.y())/height() < 0.6) {
      x = (m.x() / float(width()) * 2) - 1;
      v = 0;
    }
    return true;
  }
  bool onMouseDrag(const Mouse &m) override {
    Rayd r = getPickRay(m.x(), m.y());
    if (float(m.y())/height() > 0.4 && float(m.y())/height() < 0.6) {
      x = (m.x() / float(width()) * 2) - 1;
      v = 0;
    }
    return true;
  }
  bool onMouseUp(const Mouse &m) override {
    Rayd r = getPickRay(m.x(), m.y());
    return true;
  }
};


int main()
{
    MyApp app;
    // app.configureAudio(44100, 512, 2, 2);

    app.start();
    return 0;
}