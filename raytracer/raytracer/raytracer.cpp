// raytracer.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#define X_PIXELS 600
#define Y_PIXELS 400



#include <iostream>
#include <math.h>
#include <vector>
#include <iterator>
#include <fstream>
#include <bitset>
#include "Eigen/Dense"

using namespace std;
using namespace Eigen;

Vector4f normalize_4f(Vector4f v) {
    float norm = pow((pow(v[0], 2) + pow(v[1], 2) + pow(v[2], 2)), 0.5);
    Vector4f return_vector = { v[0] / norm, v[1] / norm, v[2] / norm, 1 };
    return return_vector;
}

Vector3f dropHomogenous(Vector4f v) {
    return Vector3f(v[0], v[1], v[2]);
}

class Color {
public:
    float red;
    float green;
    float blue;
    Color() {}
    Color(float r, float g, float b) {
        red = r;
        green = g;
        blue = b;
    }
};

class Ray {
public:
    Vector4f origin;
    Vector4f direction;
    Ray(Vector4f origin_input, Vector4f direction_input) {
        origin = origin_input;
        direction = normalize_4f(direction_input);

    }
};

class Primative {
public:
    int id;
    Color objectColor;
    Primative(Color c) {
        objectColor = c;
    }
    Primative() {}
    virtual void transform(Matrix4f m) = 0; // Transforms the 
    virtual float intersect(Ray r) = 0; // Calculates if a ray intersects the object, returns omega if it does, -1 if it does not
    virtual Vector4f position() = 0; // This should be a general indication of where the object is located. Different for each object.
};



class World {
    vector<Primative*> objectList;
    Color background = Color(0, 0, 0);

public:
    int add(Primative* obj) {

        obj->id = objectList.size();
        objectList.push_back(obj);
        return obj->id;
    }

    Color spawnRay(Ray r) {
        float closest_distance = objectList[0]->intersect(r);
        Primative* closest = objectList[0];
        for (uint16_t i = 0; i < objectList.size(); i++) {

            float intersect = objectList[i]->intersect(r);
            if (intersect > 0) {
               // cout  << intersect << ", " << closest_distance <<  endl;
            }
            if (intersect > 0  && (intersect < closest_distance || closest_distance < 0)) {

                closest = objectList[i];
                closest_distance = intersect;
            }
        }
        if (closest_distance == -1) {
            return background;
        }
        else {
            Color c = closest->objectColor;
            //im assuming more stuff will go here when shading is implemented
            return c;
        }
    }

    void transform(Matrix4f transform, int transform_id) {
        for (uint16_t i = 0; i < objectList.size(); i++) {
            if (objectList[i]->id == transform_id) {
                objectList[i]->transform(transform);
            }
        }
    }

    void transformAllObjects(Matrix4f transform) {
        for (uint16_t i = 0; i < objectList.size(); i++) {
            objectList[i]->transform(transform);
        }
    }

    Primative* getObject(int id) {
        for (uint16_t i = 0; i < objectList.size(); i++) {
            if (objectList[i]->id == id) {
                return (objectList[i]);
            }
        }
    }

};

class Camera {
public:
    Vector4f position;
    Vector4f lookat;
    Vector4f up;
    float focal_length;
    int heightPixels;
    int widthPixels;
    float planeHeight;
    float planeWidth;
    Camera(Vector4f pos, Vector4f la, Vector4f upVec, float f, int hP, int wP, float pH, float pW) {
        position = pos;
        lookat = la;
        up = upVec;
        focal_length = f;
        heightPixels = hP;
        widthPixels = wP;
        planeHeight = pH;
        planeWidth = pW;
    }
    void render(World w) {
        Matrix4f camera_transform;
        cout << "RENDERING" << endl;
        Vector3f n = (dropHomogenous(position) - dropHomogenous(lookat)).normalized();
        cout << "n is " << n << endl;
        Vector3f u = dropHomogenous(up).cross(n).normalized();
        cout << "u is " << u << endl;
        Vector3f v = n.cross(u).normalized();
        cout << "v is " << v << endl;
        Matrix4f view_trans;
        cout << "Position is " << dropHomogenous(position) << endl << endl;
        cout << dropHomogenous(position).dot(u) << endl;
        float translate_u = -dropHomogenous(position).dot(u);
        float translate_v = -dropHomogenous(position).dot(v);
        float translate_n = -dropHomogenous(position).dot(n);
        view_trans << u[0], v[0], n[0], 0,
            u[1], v[1], n[1], 0,
            u[2], v[2], n[2], 0,
            translate_u, translate_v, translate_n, 1;
        view_trans.transposeInPlace();
        cout << view_trans << endl;

        
        Matrix4f inverse_view_trans = view_trans.inverse();
        cout << inverse_view_trans * Vector4f(0, 0, focal_length, 1) << endl;
        float pixelHeight = planeHeight / heightPixels;
        float pixelWidth = planeWidth / widthPixels;
        vector<vector<Color>> image(widthPixels, vector<Color> (heightPixels));
        cout << "HERE1" << endl;
        for (uint32_t i = 0; i < widthPixels; i++) {
            for (uint32_t j = 0; j < heightPixels; j++) {
                Ray r = Ray(position, 
                    -(position - (inverse_view_trans * Vector4f((-planeWidth / 2) + (pixelWidth / 2) + (i * pixelWidth), (planeHeight / 2) - (pixelHeight / 2) - (j * pixelHeight), focal_length, 1))) );
                //cout << r.direction << endl << endl;
                Color c = w.spawnRay(r);
                image[j][i] = c;
            }
        }
        ofstream file;
        file.open("image.ppm");
        file << "P3" << endl;
        file << widthPixels << " " << heightPixels << endl;
        file << "255" << endl;

        
        for (uint32_t i = 0; i < widthPixels; i++) {
            for (uint32_t j = 0; j < heightPixels; j++) {
                file << image[i][j].red << " ";
                file << image[i][j].green<< " ";
                file << image[i][j].blue<< "\t";
            }
        }

        file.close();
    }
};

class Sphere : virtual public Primative {
public:
    Vector4f center;
    float radius;
    void transform(Matrix4f m) {
        center = m * center;
    }
    float intersect(Ray r) {
        
        float A = 1;
        float B = 2 * (r.direction[0] * (r.origin[0] - center[0]) + r.direction[1] * (r.origin[1] - center[1]) + r.direction[2] * (r.origin[2] - center[2]));
        float C = ((r.origin[0] - center[0]) * (r.origin[0] - center[0])) + 
            ((r.origin[1] - center[1]) * (r.origin[1] - center[1])) + 
            ((r.origin[2] - center[2]) * (r.origin[2] - center[2])) -
            (radius * radius);
        //cout << "Calculating intersect " << A << ", " << B << "B, " << C << ", " << (B * B) - (4 * C)<<  endl;
        if ((B * B) - (4 * C) < 0) {
            return -1;
        }
        else {
            float omega_1 = (-B + sqrt((B * B) - (4 * C))) / 2;
            float omega_2 = (-B - sqrt((B * B) - (4 * C))) / 2;
            if (omega_1 < 0 && omega_2 < 0) {
                return -1;
            }
            else if (omega_1 < 0 && omega_2 >= 0) {
                return omega_2;
            }
            else if (omega_1 >= 0 && omega_2 < 0) {
                return omega_1;
            }
            else {
                return min(omega_1, omega_2);
            }
        }
    }

    Vector4f position() {
        return center;
    }

    Sphere(Color c, Vector4f cen, float r) {
        objectColor = c;
        center = cen;
        radius = r;
    }
    

};

class Triangle : virtual public Primative {
public:
    Vector4f p0;
    Vector4f p1;
    Vector4f p2;
    void transform(Matrix4f m) {
        p0 = m * p0;
        p1 = m * p1;
        p2 = m * p2;
    }
    float intersect(Ray r) {
        Vector3f e1 = dropHomogenous(p1) - dropHomogenous(p0);
        Vector3f e2 = dropHomogenous(p2) - dropHomogenous(p0);
        Vector3f T = dropHomogenous(r.origin) - dropHomogenous(p0);
        Vector3f P = dropHomogenous(r.direction).cross(e2);
        Vector3f Q = T.cross(e1);
        if (P.dot(e1) == 0) {
            return -1;
        }
        Vector3f result_vector = (1 / (P.dot(e1))) *
            Vector3f(Q.dot(e2), P.dot(T), Q.dot(dropHomogenous(r.direction)));
        if (result_vector[0] < 0) {
            return -1;
        } else if (result_vector[1] < 0 || result_vector[2] < 0 || (result_vector[1] + result_vector[2] > 1)){
            return -1;
        } else {
            return result_vector[0];
        }
    }
    Vector4f position() {
        return p0;
    }
    Triangle(Color col, Vector4f a, Vector4f b, Vector4f c) {
        objectColor = col;
        p0 = a;
        p1 = b;
        p2 = c;
    }
};




int main()
{
    Camera* cam = new Camera(Vector4f(4.013, 2.016, 4.238, 1), Vector4f(4.013, 2.016, 3.238, 1), Vector4f(0, 1, 0, 1),
        -1, 1000, 1000, 2, 2);
    World world;
    Sphere* s = new Sphere(Color(255, 255, 0), { 3.7, 1.919, 2.58, 1 }, .5);
    Sphere* s2 = new Sphere(Color(0, 0, 255), { 4.54, 1.746, 1.41, 1 }, .5);
    
    int test_id = world.add(s);
    world.add(s2);
    //world.add(new Sphere(Color(0, 255, 0), { 5.34, 1.746, 0.41, 1 }, 1));
    world.add(new Triangle(Color(255, 0, 0),
        { 3.264 + 3.5, 1.011 , 0.052 + 5, 1 },
        { 3.264 - 1.5, 1.011, 0.052 + 5, 1 },
        { 3.264 + 3.5, 1.011, 0.052 - 5, 1 }
    ));
    world.add(new Triangle(Color(124, 0, 0),
        { 3.264 - 1.5, 1.011, 0.052 - 5, 1 },
        { 3.264 + 3.5, 1.011, 0.052 - 5, 1 },
        { 3.264 - 1.5, 1.011, 0.052 + 5, 1 }
        
    )); 
    
    cam->render(world);
}