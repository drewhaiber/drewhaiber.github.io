

#define X_PIXELS 600
#define Y_PIXELS 400
#define CHECKSIZE 0.5
#define MAX_DEPTH 10 // max ray recursion depth
#define BUMP_MAPPING 0
#define KD_DEPTH_MAX 100 //max kd tree depth
#define LIGHT_SCALE 1
#define BUNNY true
#define TONE_TYPE false //false if Ward, true if Reinhardt
#define ADV_TONE true //uses adaptive log 

#include <cmath>
#include <sstream>
#include <cstdio>
#include <ctime>
#include <fstream>
#include <chrono>
#include <iostream>
#include <math.h>
#include <vector>
#include <iterator>
#include <fstream>
#include <bitset>
#include "Eigen/Dense"
#include "quartic.h"

using namespace std;
using namespace Eigen;

//Normalizes homogenous vector based on first 3 elements
Vector4d normalize_4d(Vector4d v) {
    double norm = pow((pow(v[0], 2) + pow(v[1], 2) + pow(v[2], 2)), 0.5);
    Vector4d return_vector = { v[0] / norm, v[1] / norm, v[2] / norm, 1 };
    return return_vector;
}

//Returns the reflection of a across normal n
Vector3d reflect(Vector3d a, Vector3d n) {
    return a - (2 * (a.dot(n.normalized())) * n.normalized());
}

//converts homogenous vector to a non-homogenous vector
Vector3d dropHomogenous(Vector4d v) {
    return Vector3d(v[0], v[1], v[2]);
}

//Makes a vector homogenous
Vector4d makeHomogenous(Vector3d v) {
    return Vector4d(v[0], v[1], v[2], 1);
}



struct Voxel {
    double x_min;
    double x_max;
    double y_min;
    double y_max;
    double z_min;
    double z_max;
};

//Contains two voxels that combine into one larger voxel
struct Partitioned_Voxel {
    Voxel v1;
    Voxel v2;
};


//Returns true if a point is inside of a voxel
bool AABB_point_intersect(Vector3d point, Voxel AABB) {
    return AABB.x_min < point[0] && AABB.x_max > point[0] &&
        AABB.y_min < point[1] && AABB.y_max > point[1] &&
        AABB.z_min < point[2] && AABB.z_max > point[2];
}

//Returns true if v has any portion inside of AABB
bool AABB_intersect(Voxel v, Voxel AABB) {
    return ((v.x_min <= AABB.x_max && v.x_max >= AABB.x_min) &&
        (v.y_min <= AABB.y_max && v.y_max >= AABB.y_min) &&
        (v.z_min <= AABB.z_max && v.z_max >= AABB.z_min)) ||(
        (v.x_min >= AABB.x_min && v.x_max <= AABB.x_max) &&
        (v.y_min >= AABB.y_min && v.y_max <= AABB.y_max) &&
        (v.z_min >= AABB.z_min && v.z_max <= AABB.z_max));

}



class Color {
public:
    double red;
    double green;
    double blue;
    Color() {}
    Color(double r, double g, double b) {
        red = r;
        green = g;
        blue = b;
    }
    //Multiplies each color value by a constant
    Color mult(double s) {
        return Color(red * s, green * s, blue * s);
    }
    //Adds two colors together
    Color add(Color c) {
        return Color(red + c.red, green + c.green, blue + c.blue);
    }
    //Element wise multiplication of two colors
    Color mult_color(Color c) {
        return Color(red * c.red, green * c.green, blue * c.blue);
    }
};

struct Intersect {

    double omega; //Distance from ray origin
    Vector3d normal; //Normal at intersection
};

class Ray {
public:
    Vector3d origin;
    Vector3d direction;
    Vector3d inverses; //multiplicative inverse of each direction component for faster calculations
    double index; //Index of refraction of the medium the ray is currently in
    Ray(Vector3d origin_input, Vector3d direction_input) {
        origin = origin_input;
        direction = direction_input.normalized();
        inverses = Vector3d( 1 / direction[0], 1 / direction[1], 1 / direction[2]);
        index = 1;

    }
    Ray(Vector3d origin_input, Vector3d direction_input, double i){
        origin = origin_input;
        direction = direction_input.normalized();
        inverses = { 1 / direction[0], 1 / direction[1], 1 / direction[2] };
        index = i;
    }
};

//Abstract class, should only be instantiated as a specific primative
class Primative {
public:
    int id;
    Color objectColor;
    Color specColor;
    double ka;
    double kd;
    double ks;
    double ke;
    double kr;
    double kt;
    double index_of_refraction;

    //Gets the color of the object at the given point
    virtual Color getColor(Vector3d intersect_point) = 0;

    Primative(Color c, Color sc) {
        objectColor = c;
        specColor = sc;
    }

    Primative() {}

    virtual bool voxel_intersect(Voxel v) = 0;
    virtual void transform(Matrix4d m) = 0; // Transforms the object (NOT IMPLEMENTED, AS IT WAS UNNECESSARY)
    virtual Intersect intersect(Ray r) = 0; // Calculates if a ray intersects the object, returns omega if it does, -1 if it does not
    virtual Vector3d position() = 0; // This should be a general indication of where the object is located. Different for each object.
};


//KD tree nodes
class Node {
public:
    bool is_leaf;
    Node* front;
    Node* rear;
    vector<Primative*> objects;
    int axis; //0 = x, 1 = y, 2 = z
    double plane;
    Voxel box;

    //Calculates the size of the tree for debugging
    int size() {
        if (is_leaf) {
            return 1;
        }
        else {
            return 1 + front->size() + rear->size();
        }
    }
    //Non-leaf node constructor
    Node(vector<Primative*> objects_list, Node* front_node, Node* rear_node, int axis_type, double plane_value, Voxel voxel) {
        front = front_node;
        rear = rear_node;
        axis = axis_type;
        plane = plane_value;
        is_leaf = false;
        box = voxel;
        objects = objects_list;
    }


    //Leaf node constructor
    Node(vector<Primative*> objects_list, Voxel voxel) {
        objects = objects_list;
        is_leaf = true;
        box = voxel;
    }

    //Deletes the node and sub-nodes
    void delete_rec() {
        if (!is_leaf) {
            front->delete_rec();
            rear->delete_rec();
        }
        free(this);
    }
    
};



//Partitions a voxel into sub-voxels across axis "type" at "value"
Partitioned_Voxel partition(Voxel v, int type, double value) {
    Voxel v1, v2;
    if (type == 0){
        v1 = { value, v.x_max, v.y_min, v.y_max, v.z_min, v.z_max };
        v2 = { v.x_min, value, v.y_min, v.y_max, v.z_min, v.z_max };
    }
    else if (type == 1) {
        v1 = { v.x_min, v.x_max, value, v.y_max, v.z_min, v.z_max };
        v2 = { v.x_min, v.x_max, v.y_min, value, v.z_min, v.z_max };
    }
    else {
        v1 = { v.x_min, v.x_max, v.y_min, v.y_max, value, v.z_max };
        v2 = { v.x_min, v.x_max, v.y_min, v.y_max, v.z_min, value };
    }
    Partitioned_Voxel to_return = { v1, v2 };
    return to_return;
}

//Determines when KD tree stops generating deeper nodes
bool terminal(vector<Primative*> objects, Voxel v, int depth) {
    return (objects.size() <= 2 || (v.x_max - v.x_min < .0005 || v.y_max - v.y_min < .0005 || v.z_max - v.z_min < .0005)  ||
         depth > KD_DEPTH_MAX);
}


//Constructs the tree
Node* getNode(vector<Primative*> objects, Voxel v, int prev_axis, int depth) {
    //Endpoint reached
    if (terminal(objects, v, depth)) {
        vector<Primative*> new_objects = objects;
        return new Node(new_objects, v);
    }
    
    //Determines split point
    double p = 0;
    if (prev_axis == 0) {
        p = (v.x_min + v.x_max) * .5;
    }
    else if (prev_axis == 1) {
        p = (v.y_min + v.y_max) * .5;
    }
    else {
        p = (v.z_min + v.z_max) * .5;
    }
    
    //Partition voxel into sub-voxels
    Partitioned_Voxel partitioned = partition(v, prev_axis, p);
    Voxel v1 = partitioned.v1;
    Voxel v2 = partitioned.v2;
    vector<Primative*> objects1;
    vector<Primative*> objects2;
    for (int i = 0; i < objects.size(); i++) {
        if (objects[i]->voxel_intersect(v1)) {
            objects1.push_back(objects[i]);
        }
        if (objects[i]->voxel_intersect(v2)) {
            objects2.push_back(objects[i]);
        }
    }
    //Get sub Nodes
    return new Node(objects, getNode(objects1, v1, (prev_axis + 1)%3, depth+1), getNode(objects2, v2, (prev_axis + 1)%3, depth+1), prev_axis, p, v);
}


struct IntersectInfo {
    Intersect intersect;
    Primative* object;
};

class World {

public:

    vector<Primative*> objectList;
    Color background = Color(.1 *LIGHT_SCALE, .1*LIGHT_SCALE , .5*LIGHT_SCALE );
    Vector3d light_position = { 4, 5, 3.7 };
    Color light_color = Color(LIGHT_SCALE, LIGHT_SCALE, LIGHT_SCALE);
    Node* tree;
    //Adds a primative to the scene
    int add(Primative* obj) {

        obj->id = objectList.size();
        objectList.push_back(obj);
        return obj->id;
    }

    //Determines if a voxel and ray intersect, if it does the distance from origin
    //is at tmin for first intersection, tmax for second intersection
    bool intersection(Voxel v, Ray r, double* tmin, double* tmax) {
        double tx1 = (v.x_min - r.origin[0]) * r.inverses[0];
        double tx2 = (v.x_max - r.origin[0]) * r.inverses[0];
        *tmin = min(tx1, tx2);
        *tmax = max(tx1, tx2);

        double ty1 = (v.y_min - r.origin[1]) * r.inverses[1];
        double ty2 = (v.y_max - r.origin[1]) * r.inverses[1];

        *tmin = max(*tmin, min(ty1, ty2));
        *tmax = min(*tmax, max(ty1, ty2));

        double tz1 = (v.z_min - r.origin[2]) * r.inverses[2];
        double tz2 = (v.z_max - r.origin[2]) * r.inverses[2];

        *tmin = max(*tmin, min(tz1, tz2));
        *tmax = min(*tmax, max(tz1, tz2));

        return *tmax >= *tmin;
    }

    //Finds the intersection of a ray with primatives in the Kd tree "node"
    IntersectInfo traverse(Ray r, Node* node, int axis) {
        //Tests intersect against elements within leaf node
        if (node->is_leaf) {
            if (node->objects.size() == 0){
                IntersectInfo to_return;
                to_return.intersect = { -1, {0, 0, 0} };
                to_return.object = NULL;
                return to_return;
            }
            Intersect closest_intersect = node->objects[0]->intersect(r);
            double closest_distance = closest_intersect.omega;
            Primative* closest = node->objects[0];
            for (uint16_t i = 0; i < node->objects.size(); i++) {

                Intersect current_intersect = node->objects[i]->intersect(r);
                double intersect = current_intersect.omega;
                if (intersect > 0 && (intersect < closest_distance || closest_distance < 0)) {
                    closest = node->objects[i];
                    closest_distance = intersect;
                    closest_intersect = current_intersect;
                }
            }
            IntersectInfo to_return;
 
            to_return.intersect = closest_intersect;
            to_return.object = closest;
            return to_return;
        }

        //Traverses node recursively
        Vector3d a, b;
        double t_min;
        double t_max;
        intersection(node->box, r, &t_min, &t_max);

        if (t_min < 0) { //ray is inside the box
            a = r.origin;
        }
        else {
            a = r.origin + (r.direction * t_min);
        }
        b = r.origin + (r.direction * t_max);
        
        //Determine which order of subnodes to traverse
        if (a[axis] <= node->plane) {
            if (b[axis] < node->plane) {
                return traverse(r, node->rear, (axis + 1) % 3);
            }
            else {
                IntersectInfo left = traverse(r, node->rear, (axis + 1) % 3);
                if (left.intersect.omega < 0) {
                    return traverse(r, node->front, (axis + 1) % 3);
                }
                else {
                    return left;
                }
            }
        }
        else {
            if (b[axis] <= node->plane) {

                IntersectInfo right = traverse(r, node->front, (axis + 1) % 3);
                if (right.intersect.omega < 0) {
                    return traverse(r, node->rear, (axis + 1) % 3);
                }
                else {
                    return right;
                }
            }
            else {
                return traverse(r, node->front, (axis + 1) % 3);
            }
        }
    }
    //Spawns a ray into the world
    Color spawnRay(Ray r, int depth) {
        IntersectInfo info = traverse(r, tree, tree->axis);
        //FIND CLOSEST OBJECT
        Intersect closest_intersect = info.intersect;
        double closest_distance = closest_intersect.omega;
        Primative* closest = info.object;

        if (closest_distance < 0) {
            //No intersection
            return background;
        }
        else {
            //Intersection
            Color objectColor = closest->getColor(r.origin + (r.direction * closest_intersect.omega));
            Color c = objectColor.mult(closest->ka); //ambient color 
            Vector3d shadow_direction = (Vector3d(light_position[0], light_position[1], light_position[2]) - (r.origin + (closest_distance * r.direction))).normalized();
            Ray shadowRay = Ray(r.origin + (closest_distance * r.direction) + (.001 * Vector3d(shadow_direction[0], shadow_direction[1], shadow_direction[2])), shadow_direction);
            IntersectInfo shadow_info = traverse(shadowRay, tree, tree->axis);
            uint16_t sees_light = 1;
            bool has_transparent = false;
            while (shadow_info.intersect.omega > 0) {
                //Detects if there is a shadow or not
                if (shadow_info.object->kt != 0) {
                    has_transparent = true;
                    Vector3d new_shadow_direction = (light_position - shadowRay.origin).normalized();
                    shadowRay = Ray(shadowRay.origin + (shadow_info.intersect.omega * shadow_direction) + (.001 * shadow_direction), new_shadow_direction);
                    shadow_info = traverse(shadowRay, tree, tree->axis);
                }
                else {
                    
                    sees_light = 0;
                    break;
                }
            }

            Vector3d reflect_vector = reflect(shadowRay.direction, closest_intersect.normal).normalized();
            if(sees_light == 1){
                if (has_transparent) {
                    c = c.add((objectColor).mult_color(light_color).mult(max(0.0, shadowRay.direction.dot(closest_intersect.normal))).mult(closest->kd/2));

                    c = c.add((closest->specColor).mult_color(light_color).mult(.5 * closest->ks * pow(max(0.0, reflect_vector.dot(r.direction)), closest->ke)));
                }
                else {
                    c = c.add((objectColor).mult_color(light_color).mult(max(0.0, shadowRay.direction.dot(closest_intersect.normal))).mult(closest->kd ));

                    c = c.add((closest->specColor).mult_color(light_color).mult( closest->ks * pow(max(0.0, reflect_vector.dot(r.direction)), closest->ke)));
                }
                
            }

            if (depth < MAX_DEPTH) {
                if (closest->kr > 0) {// If object is reflective
                    Vector3d reflect_rec = reflect(r.direction, closest_intersect.normal);
                    Ray r_reflect = Ray(r.origin + (r.direction * closest_intersect.omega) + (.001*reflect_rec), reflect_rec);
                    c = c.add(this->spawnRay(r_reflect, depth + 1).mult(closest->kr));
                }
                if (closest->kt > 0) { //If object is transparent
                    Vector3d normal = closest_intersect.normal.normalized();
                    double cos_i = normal.dot(-r.direction.normalized());
                    double index_ratio = r.index / closest->index_of_refraction;
                    if (cos_i < 0) { //If backside of object
                        normal = -normal;
                        cos_i = normal.dot(-r.direction.normalized());
                        index_ratio = 1 / r.index;
                    }
                    double sin2t = pow(index_ratio, 2) * (1 - pow(cos_i, 2));
                    double rad = 1 - sin2t;
                    if (rad < 0) {//Total internal reflection
                        Vector3d reflect_rec = reflect(r.direction, normal);
                        Ray r_reflect = Ray(r.origin + (r.direction * closest_intersect.omega) + (.001 * reflect_rec), reflect_rec);
                        c = c.add(this->spawnRay(r_reflect, depth + 1).mult(closest->kt));
                    } else { //Not internal reflection
                        Vector3d trans_rec = ((index_ratio) * r.direction.normalized()) +
                            (((r.index / closest->index_of_refraction) *cos_i) - sqrt(rad)) * (normal);
                        Ray r_trans = Ray(r.origin + (r.direction * closest_intersect.omega) + (.001 * trans_rec), trans_rec, 1);
                        if ((closest->index_of_refraction - r.index) < 0.00001 && (closest->index_of_refraction - r.index) > -0.00001) {
                            //Going from object to air
                            r_trans = Ray(r.origin + (r.direction * closest_intersect.omega) + (.001 * trans_rec), trans_rec, 1);
                        } else {
                            //Going from air to object
                            r_trans = Ray(r.origin + (r.direction * closest_intersect.omega) + (.001 * trans_rec), trans_rec, closest->index_of_refraction);
                          }
                        Color temp = this->spawnRay(r_trans, depth + 1).mult(closest->kt);
                        if (cos_i < 0) {
                            temp = temp.mult(1 / closest->kt);
                        }
                        c = c.add(temp);
                    }
                    
                }
            }
            return c;
        }
    }

    //Unused
    void transform(Matrix4d transform, int transform_id) {
        for (uint16_t i = 0; i < objectList.size(); i++) {
            if (objectList[i]->id == transform_id) {
                objectList[i]->transform(transform);
            }
        }
    }

    //Unused
    void transformAllObjects(Matrix4d transform) {
        for (uint16_t i = 0; i < objectList.size(); i++) {
            objectList[i]->transform(transform);
        }
    }

    //Gets the object with id
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
    Vector3d position;
    Vector3d lookat;
    Vector3d up;
    double focal_length;
    int heightPixels;
    int widthPixels;
    double planeHeight;
    double planeWidth;
    Camera(Vector3d pos, Vector3d la, Vector3d upVec, double f, int hP, int wP, double pH, double pW) {
        position = pos;
        lookat = la;
        up = upVec;
        focal_length = f;
        heightPixels = hP;
        widthPixels = wP;
        planeHeight = pH;
        planeWidth = pW;
    }
    //Renders the scene
    void render(World w) {
        Matrix4d camera_transform;
        cout << "RENDERING" << endl;
        Vector3d n = (position - lookat).normalized();
        Vector3d u = up.cross(n).normalized();
        Vector3d v = n.cross(u).normalized();
        Matrix4d view_trans;

        double translate_u = -position.dot(u);
        double translate_v = -position.dot(v);
        double translate_n = -position.dot(n);
        view_trans << u[0], v[0], n[0], 0,
            u[1], v[1], n[1], 0,
            u[2], v[2], n[2], 0,
            translate_u, translate_v, translate_n, 1;
        view_trans.transposeInPlace();


        Matrix4d inverse_view_trans = view_trans.inverse();
        double pixelHeight = planeHeight / heightPixels;
        double pixelWidth = planeWidth / widthPixels;
        vector<vector<Color>> image(widthPixels, vector<Color>(heightPixels));
        Vector4d position4 = makeHomogenous(position);
        for (int i = 0; i < widthPixels; i++) {
            for (int j = 0; j < heightPixels; j++) {

                //4 rays for supersampling
                Ray r1 = Ray(position,
                    dropHomogenous(-(position4 - (inverse_view_trans * Vector4d((-planeWidth / 2) + (pixelWidth / 4) + (i * pixelWidth), (planeHeight / 2) - (pixelHeight / 4) - (j * pixelHeight), focal_length, 1)))));
                Ray r2 = Ray(position,
                    dropHomogenous(-(position4 - (inverse_view_trans * Vector4d((-planeWidth / 2) + (3*pixelWidth / 4) + (i * pixelWidth), (planeHeight / 2) - (pixelHeight / 4) - (j * pixelHeight), focal_length, 1)))));
                Ray r3 = Ray(position,
                    dropHomogenous(-(position4 - (inverse_view_trans * Vector4d((-planeWidth / 2) + (3*pixelWidth / 4) + (i * pixelWidth), (planeHeight / 2) - (3*pixelHeight / 4) - (j * pixelHeight), focal_length, 1)))));
                Ray r4 = Ray(position,
                    dropHomogenous(-(position4 - (inverse_view_trans * Vector4d((-planeWidth / 2) + (pixelWidth / 4) + (i * pixelWidth), (planeHeight / 2) - (3*pixelHeight / 4) - (j * pixelHeight), focal_length, 1)))));

                Color c = w.spawnRay(r1, 0).add(w.spawnRay(r2, 0)).add(w.spawnRay(r3, 0)).add(w.spawnRay(r4, 0)).mult(.25);
                
                image[j][i] = c;
            }
        }

        vector<vector<double>> L(widthPixels, vector<double>(heightPixels));
        double max = 0.f;
        double log_sum = 0.f;
        double epsilon = 1.f / 1000000.f;
        float ldmax = 1;

        for (uint32_t i = 0; i < widthPixels; i++) {
            for (uint32_t j = 0; j < heightPixels; j++){
                L[i][j] = (image[i][j].red * .27) + (image[i][j].green * .67) + (image[i][j].blue * .06);
                max = max > L[i][j] ? max : L[i][j];
                log_sum += log(epsilon + L[i][j]);
            }
        }

        double log_luminance = exp((1.f / (widthPixels * heightPixels)) * (log_sum));

        double sf = pow((1.219 + pow(((float)max)/2, .4))/(1.219 + pow(log_luminance, .4)), 2.5);

        ofstream file;
        file.open("image.ppm");
        file << "P3" << endl;
        file << widthPixels << " " << heightPixels << endl;
        int scale = 255;//Numbers are normally 0-LIGHT_SCALE, but ppm only accepts integers
                        //Scale makes the new scale integers between 0 and LIGHT_SCALE * scale for each color.
                        //Useful for low lighting environments

        file << round(ldmax * scale) << endl;

        float bias = .85;
        if (ADV_TONE) {
            for (uint32_t i = 0; i < widthPixels; i++) {
                for (uint32_t j = 0; j < heightPixels; j++) {
                    float denom = log(2 + (pow((L[i][j]/log_luminance) / (max/log_luminance), (log(bias) / log(0.5)))*8));
                    float new_luminance = (1 / log10(1 + (max / log_luminance))) * log(1 + (L[i][j] / log_luminance)) / denom;
                    float scale_factor = new_luminance / L[i][j];
                    file << round(scale * min(scale_factor * image[i][j].red, (double)ldmax)) << " ";
                    file << round(scale * min(scale_factor * image[i][j].green, (double)ldmax)) << " ";
                    file << round(scale * min(scale_factor * image[i][j].blue, (double)ldmax)) << " ";

                }
            }
        } else {
            if (!TONE_TYPE) {//Ward tone
                for (uint32_t i = 0; i < widthPixels; i++) {
                    for (uint32_t j = 0; j < heightPixels; j++) {
                        file << round(scale * min(sf * image[i][j].red, (double)ldmax)) << " ";
                        file << round(scale * min(sf * image[i][j].green, (double)ldmax)) << " ";
                        file << round(scale * min(sf * image[i][j].blue, (double)ldmax)) << "\t";
                    }
                }
            }
            else {//Reinhard tone
                vector<vector<Color>> scaled(widthPixels, vector<Color>(heightPixels));
                vector<vector<Color>> reflectance(widthPixels, vector<Color>(heightPixels));

                double a = 0.18;
                double key_value = log_luminance;
                for (uint32_t i = 0; i < widthPixels; i++) {
                    for (uint32_t j = 0; j < heightPixels; j++) {
                        scaled[i][j] = image[i][j].mult(a / log_luminance);
                        reflectance[i][j] = Color(scaled[i][j].red / (1 + scaled[i][j].red),
                            scaled[i][j].green / (1 + scaled[i][j].green),
                            scaled[i][j].blue / (1 + scaled[i][j].blue));
                        file << round(scale * min(reflectance[i][j].red * ldmax, (double)ldmax)) << " ";
                        file << round(scale * min(reflectance[i][j].green * ldmax, (double)ldmax)) << " ";
                        file << round(scale * min(reflectance[i][j].blue * ldmax, (double)ldmax)) << "\t";
                    }
                }

            }
        }

        file.close();
    }
};

class Sphere : virtual public Primative {
public:
    Vector3d center;
    double radius;

    Color getColor(Vector3d intersect_point) {
        return objectColor;
    }

    bool voxel_intersect(Voxel v) {

        double d = 0;
        
        if (center[0] < v.x_min) {
            d += (center[0] - v.x_min)*(center[0] - v.x_min);
        }
        else if (center[0] > v.x_max) {
            d += (center[0] - v.x_max) * (center[0] - v.x_max);
        }
        if (center[1] < v.y_min) {
            d += (center[1] - v.y_min) * (center[1] - v.y_min);
        }
        else if (center[1] > v.y_max) {
            d += (center[1] - v.y_max) * (center[1] - v.y_max);
        }
        if (center[2] < v.z_min) {
            d += (center[2] - v.z_min) * (center[2] - v.z_min);
        }
        else if (center[2] > v.z_max) {
            d += (center[2] - v.z_max) * (center[2] - v.z_max);
        }
        return d <= radius * radius;
    }

    void transform(Matrix4d m) {
        center = dropHomogenous(m * makeHomogenous(center));
    }
    Intersect intersect(Ray r) {
        
        double A = 1;
        double B = 2 * (r.direction[0] * (r.origin[0] - center[0]) + r.direction[1] * (r.origin[1] - center[1]) + r.direction[2] * (r.origin[2] - center[2]));
        double C = ((r.origin[0] - center[0]) * (r.origin[0] - center[0])) + 
            ((r.origin[1] - center[1]) * (r.origin[1] - center[1])) + 
            ((r.origin[2] - center[2]) * (r.origin[2] - center[2])) -
            (radius * radius);
        if ((B * B) - (4 * C) < 0) {
            return { -1, Vector3d(0, 0, 0) };
        }
        else {

            double omega_1 = (-B + sqrt((B * B) - (4 * C))) / 2;
            double omega_2 = (-B - sqrt((B * B) - (4 * C))) / 2;
            if (omega_1 < 0 && omega_2 < 0) {
                return { -1, Vector3d(0, 0, 0) };
            }
            else if (omega_1 < 0 && omega_2 >= 0) {
                return { omega_2, ((r.origin + (omega_2 * r.direction)) - center).normalized() };

            }
            else if (omega_1 >= 0 && omega_2 < 0) {
                return { omega_1, ((r.origin + (omega_1 * r.direction)) - center).normalized() };
            }
            else {
                return omega_1 < omega_2 ?
                    Intersect{ omega_1, ((r.origin + (omega_1 * r.direction))
                        - center).normalized() } :
                Intersect {omega_2, ((r.origin + (omega_2 * r.direction))
                    - center).normalized()};
            }
        }
    }

    Vector3d position() {
        return center;
    }

    Sphere(Color c, Color sc, Vector3d cen, double r, double KE, double KR, double KD, double KS, double KT, double index) {
        objectColor = c;
        specColor = sc;
        center = cen;
        radius = r;
        ka = .1 * LIGHT_SCALE;
        kd = KD;
        ks = KS;
        ke = KE;
        kr = KR;
        kt = KT;
        index_of_refraction = index;
    }
    

};

class Torus : virtual public Primative {
public:
    Vector3d center;
    double minor_radius;
    double major_radius;

    Color getColor(Vector3d intersect_point) {
        return objectColor;
    }


    bool voxel_intersect(Voxel v) {

        Voxel box;
        box.x_min = center[0] - (major_radius + minor_radius);
        box.x_max = center[0] + (major_radius + minor_radius);
        box.y_min = center[1] - (minor_radius);
        box.y_max = center[1] + (minor_radius);
        box.z_min = center[2] - (major_radius + minor_radius);
        box.z_max = center[2] + (major_radius + minor_radius);
        return AABB_intersect(box, v);

    }

    Torus(Color c, Color sc,  Vector3d cen, double r, double R, double KR = 0, double KT = 0) {
        objectColor = c;
        specColor = sc;
        center = cen;
        minor_radius = r;
        major_radius = R;
        ka = .1 * LIGHT_SCALE;
        kd = .1;
        ks = .1;
        ke = 20;
        kr = KR;
        kt = .7;
        index_of_refraction = 1.15;
    }

    Intersect intersect(Ray r) {
        r.origin = r.origin - center;
        double c4 = pow(pow(r.direction[0], 2) + pow(r.direction[1], 2) + pow(r.direction[2], 2), 2);

        double c3 = 4 * (pow(r.direction[0], 2) + pow(r.direction[1], 2) + pow(r.direction[2], 2)) * 
            (((double)r.direction[0] * (double)r.origin[0]) + ((double)r.direction[1] * (double)r.origin[1]) + ((double)r.direction[2] * (double)r.origin[2]));

        double c2 = 2 * (pow(r.direction[0], 2) + pow(r.direction[1], 2) + pow(r.direction[2], 2)) *
            (pow(r.origin[0], 2) + pow(r.origin[1], 2) + pow(r.origin[2], 2) - pow(minor_radius, 2) - pow(major_radius, 2)) +
            (4 * pow((((double)r.direction[0] * (double)r.origin[0]) + ((double)r.direction[1] * (double)r.origin[1]) + ((double)r.direction[2] * (double)r.origin[2])), 2)) + (4 * pow(major_radius, 2) * pow(r.direction[1], 2));

        double c1 = 4 * ((((double)r.direction[0] * (double)r.origin[0]) + ((double)r.direction[1] * (double)r.origin[1]) + ((double)r.direction[2] * (double)r.origin[2])))
            * (pow(r.origin[0], 2) + pow(r.origin[1], 2) + pow(r.origin[2], 2) - pow(minor_radius, 2) - pow(major_radius, 2)) + 8 *(pow(major_radius, 2) * r.origin[1] * r.direction[1]);

        double c0 = pow((pow(r.origin[0], 2) + pow(r.origin[1], 2) + pow(r.origin[2], 2) - pow(minor_radius, 2) - pow(major_radius, 2)), 2) - (4 * pow(major_radius, 2) * (pow(minor_radius, 2) - pow(r.origin[1], 2)));

        DComplex* roots = solve_quartic(c3 / c4, c2 / c4, c1 / c4, c0 / c4);
        double min_distance = 1.7976931348623157E+308; // biggest double i believe... probably could be better implemented
        for (int i = 0; i < 4; i++) {
            if (roots[i].imag() == 0) {
                if (roots[i].real() < min_distance && roots[i].real() >= 0) {
                    min_distance = roots[i].real();
                }
            }
        }
        if (min_distance > 1.79E+307) {// this is really janky and should be replaced ASAP
            return Intersect{-1, Vector3d(1, 0, 0)};
        }
        else {
            
            Vector3d point = (r.origin + (min_distance * r.direction.normalized()));
            double alpha = major_radius / sqrt(pow(point[0], 2) + pow(point[2], 2));
            return Intersect{ min_distance, Vector3d((1 - alpha) * point[0],  point[1], (1 - alpha) * point[2]).normalized() };
        }
    }

    void transform(Matrix4d m) {
        center = dropHomogenous(m * makeHomogenous(center));
    }

    Vector3d position() {
        return center;
    }
};
class Triangle : virtual public Primative {
public:
    Vector3d p0;
    Vector3d p1;
    Vector3d p2;
    bool tex;

    bool voxel_intersect(Voxel v) {
        Voxel box;

        //Gotten from https://stackoverflow.com/questions/17458562/efficient-aabb-triangle-intersection-in-c-sharp
        //Calculates if there is a separating plane between the triangle and the voxel
        Vector3d center = Vector3d((v.x_max + v.x_min) * .5, (v.y_max + v.y_min) * .5, (v.z_max + v.z_min) * .5);
        Vector3d extents = Vector3d((v.x_max - v.x_min) * .5, (v.y_max - v.y_min) * .5, (v.z_max - v.z_min) * .5);

        Vector3d point0 = p0 - center;
        Vector3d point1 = p1 - center;
        Vector3d point2 = p2 - center;


        //edge ras
        Vector3d f[3]{
            point1 - point0,
            point2 - point1,
            point0 - point2
        };

        Vector3d box_normals[3] = {
            Vector3d(1, 0, 0),
            Vector3d(0, 1, 0),
            Vector3d(0, 0, 1)
        };
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                Vector3d axis = box_normals[i].cross(f[j]);
                float proj0 = point0.dot(axis);
                float proj1 = point1.dot(axis);
                float proj2 = point2.dot(axis);
                float r = (extents[0] * abs(box_normals[0].dot(axis))) +
                    (extents[1] * abs(box_normals[1].dot(axis))) +
                    (extents[2] * abs(box_normals[2].dot(axis)));
                if (max(-max(max(proj0, proj1), proj2), min(min(proj0, proj1), proj2)) > r) {
                    return false;
                }
            }
        }
        for (int i = 0; i < 3; i++) {
            Vector3d axis = box_normals[0];
            float proj0 = point0.dot(axis);
            float proj1 = point1.dot(axis);
            float proj2 = point2.dot(axis);
            float r = (extents[0] * abs(box_normals[0].dot(axis))) +
                (extents[1] * abs(box_normals[1].dot(axis))) +
                (extents[2] * abs(box_normals[2].dot(axis)));
            if (max(-max(max(proj0, proj1), proj2), min(min(proj0, proj1), proj2)) > r) {
                return false;
            }
        }
        Vector3d axis = f[0].cross(f[1]);
        float proj0 = point0.dot(axis);
        float proj1 = point1.dot(axis);
        float proj2 = point2.dot(axis);
        float r = (extents[0] * abs(box_normals[0].dot(axis))) +
            (extents[1] * abs(box_normals[1].dot(axis))) +
            (extents[2] * abs(box_normals[2].dot(axis)));
        if (max(-max(max(proj0, proj1), proj2), min(min(proj0, proj1), proj2)) > r) {
            return false;
        }
        return true;
    }

    Color getColor(Vector3d intersect_point) {
        if (!tex) {
            return objectColor;
        }
        double x_origin = 3.264 - 1.5;
        double z_origin = 0.052 - 5;
        double x_point = intersect_point[0];
        double z_point = intersect_point[2];
        bool col = ((int)((x_point - x_origin) / CHECKSIZE)) % 2 == 0;
        bool row = ((int)((z_point - z_origin) / CHECKSIZE)) % 2 == 0;
        if (row == col) {
            return Color(1, 0, 0);
        }
        else {
            return Color(0, 1, 0);
        }
    }

    void transform(Matrix4d m) {
        p0 = dropHomogenous(m * makeHomogenous(p0));
        p1 = dropHomogenous(m * makeHomogenous(p1));
        p2 = dropHomogenous(m * makeHomogenous(p2));
    }
    Intersect intersect(Ray r) {
        Vector3d e1 = p1 - p0;
        Vector3d e2 = p2 - p0;
        Vector3d T = r.origin - p0;
        Vector3d P = r.direction.cross(e2);
        Vector3d Q = T.cross(e1);
        if (P.dot(e1) == 0) {
            return Intersect{ -1, Vector3d(0, 0, 0) };
        }
        Vector3d result_vector = (1 / (P.dot(e1))) *
            Vector3d(Q.dot(e2), P.dot(T), Q.dot(r.direction));
        if (result_vector[0] < 0) {
            return Intersect{ -1, Vector3d(0, 0, 0) };
        } else if (result_vector[1] < 0 || result_vector[2] < 0 || (result_vector[1] + result_vector[2] > 1)){
            return Intersect{ -1, Vector3d(0, 0, 0) };
        } else {
            Vector3d normal = (e2.cross(e1)).normalized();
            Vector3d point = r.origin + (result_vector[0] * r.direction);
            double x_origin = 3.264 - 1.5;
            double z_origin = 0.052 - 5;
            double x_point = point[0];
            double z_point = point[2];
            double x_inner_coord = (x_point - x_origin) / CHECKSIZE;
            double z_inner_coord = (z_point - z_origin) / CHECKSIZE;
            double x_dif = x_inner_coord - (int)x_inner_coord;
            double z_dif = z_inner_coord - (int)z_inner_coord;
            double depth = 1; //How deep the bumps should be
            if (BUMP_MAPPING) {//Adjusts the normal if bump mapping
                if (x_dif < .15) {
                    if (z_dif < .15) {
                        if (x_dif < z_dif) {
                            normal = normal + Vector3d(depth, 0, 0);
                        }
                        else {
                            normal = normal + Vector3d(0, 0, depth);
                        }
                    }
                    else if (z_dif > .85) {
                        if (x_dif < (1 - z_dif)) {
                            normal = normal + Vector3d(depth, 0, 0);
                        }
                        else {
                            normal = normal - Vector3d(0, 0, depth);
                        }
                    }
                    else {
                        normal = normal + Vector3d(depth, 0, 0);
                    }

                }
                else if (x_dif > .85) {
                    if (z_dif < .15) {
                        if ((1 - x_dif) < z_dif) {
                            normal = normal - Vector3d(depth, 0, 0);
                        }
                        else {
                            normal = normal + Vector3d(0, 0, depth);
                        }
                    }
                    else if (z_dif > .85) {
                        if ((1 - x_dif) < (1 - z_dif)) {
                            normal = normal - Vector3d(depth, 0, 0);
                        }
                        else {
                            normal = normal - Vector3d(0, 0, depth);
                        }
                    }
                    else {
                        normal = normal - Vector3d(depth, 0, 0);
                    }

                }
                else if (z_dif < .15) {
                    normal = normal + Vector3d(0, 0, depth);
                }
                else if (z_dif > .85) {
                    normal = normal - Vector3d(0, 0, depth);
                }
            }
            
            return Intersect{ result_vector[0], normal.normalized() };
        }
    }
    Vector3d position() {
        return p0;
    }
    Triangle(Color col, Color sc, Vector3d a, Vector3d b, Vector3d c, double KR, double KT, bool texture = false) {
        objectColor = col;
        specColor = sc;
        p0 = a;
        p1 = b;
        p2 = c;
        ka = .1 * LIGHT_SCALE;
        kd = .6;
        ks = .3;
        ke = 1;
        kr = KR;
        kt = KT;
        index_of_refraction = 1.4;
        tex = texture;
    }
};

//Gets a number from a string
int get_num(string s) {
    for (int i = 0; i < s.length(); i++) {
        if (isdigit(s[i])) {
            return atoi(s.substr(i, s.length() - 1).c_str());
        }
    }
}

int main(int argc, char** argv) 
{

    Camera* cam = new Camera(Vector3d(4.013, 2.016, 4.238), Vector3d(4.013, 2.016, 3.238), Vector3d(0, 1, 0),
        -1, 1080, 1080, 2, 2);
    World world;
    
    world.add(new Sphere(Color(0, 0, .2), Color(1, 1, 1), { 3.7, 1.919, 2.58 }, .5, 50, 0, .05, .1, .75, .98));
    world.add(new Sphere(Color(1, 1, 1), Color(1, 1, 1), { 4.54, 1.746, 1.41 }, .5, 100, .75, .05, .1, 0, 1));
    world.add(new Torus(Color(0, 0.5, 1), Color(1, 1, 1), { 4.7, 1.5, 3.08 }, .1, .45, 0));
    world.add(new Triangle(Color(.5, 0, 0), Color(1, 1, 1),
        { 3.264 + 3.5, 1.011 , 0.052 + 5 },
        { 3.264 - 1.5, 1.011, 0.052 + 5 },
        { 3.264 + 3.5, 1.011, 0.052 - 5 },
        0, 0, true
    ));
    world.add(new Triangle(Color(0, .5, 0), Color (1, 1, 1),
        { 3.264 - 1.5, 1.011, 0.052 - 5 },
        { 3.264 + 3.5, 1.011, 0.052 - 5 },
        { 3.264 - 1.5, 1.011, 0.052 + 5 },
        0, 0, true
    ));

    if (BUNNY) {
        //Reads in bunny from file
        ifstream stream("bun_zipper.ply");
        string line;
        Vector3d location = {3.4, 1.2, 3.15};
        double size_scale = 3;
        if (stream.is_open()) {
            cout << "OPENED BUNNY FILE" << endl;
            int vertex_num = 0;
            int face_num = 0;
            while (getline(stream, line))
            {   
                if (line.find("element vertex") != string::npos){
                    vertex_num = get_num(line);
                }
                else if (line.find("element face") != string::npos) {
                    face_num = get_num(line);
                }
                else if (line.find("end_header") != string::npos) {
                    break;
                }
            }
            std::size_t st;
            vector<Vector3d> points(vertex_num);
            for (int i = 0; i < vertex_num; i++) {
                getline(stream, line);
                int ind = line.find(" ");
                string token = line.substr(0, ind);
                double x = stod(token);
                line = line.substr(ind+1);
                ind = line.find(" ");
                token = line.substr(0, ind);


                double y = stod(token);
                line = line.substr(ind+1);
                ind = line.find(" ");
                token = line.substr(0, ind);

                double z = stod(token);

                
                points[i] = { x, y, z };
            }
            for (int i = 0; i < face_num; i++) {
                getline(stream, line);
                int ind = line.find(" ");
                line = line.substr(ind + 1);
                ind = line.find(" ");
                string token = line.substr(0, ind);


                int p1 = stoi(token);
                line = line.substr(ind + 1);
                ind = line.find(" ");
                token = line.substr(0, ind);
                int p2 = stoi(token);
                line = line.substr(ind + 1);
                ind = line.find(" ");
                token = line.substr(0, ind);
                int p3 = stoi(token);
                line = line.substr(ind + 1);
                ind = line.find(" ");
                token = line.substr(0, ind);

                
                world.add(new Triangle(Color(0, 1, 1), Color(1, 1, 1), (size_scale * points[p1]) + location, (size_scale * points[p2]) + location,(size_scale * points[p3]) + location, 0, 0));
            }
            
            stream.close();
            cout << "CLOSE BUNNY FILE" << endl;
        }
    }

    Voxel box;
    box.x_min = 0;
    box.x_max = 8;
    box.y_min = 1;
    box.y_max = 4;
    box.z_min = -5;
    box.z_max = 7;
    cout << "STARTING TREE BUILD " << endl;
    world.tree = getNode(world.objectList, box, 0, 1);
    cout << "ENDING TREE BUILD " << endl;

    auto start = std::chrono::high_resolution_clock::now();
    cam->render(world);
    cout << "DONE RENDERING";
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::microseconds>(stop - start);
    cout << "Time taken by function: "
        << duration.count() << " microseconds" << endl;
    

    world.tree->delete_rec();


    delete cam;
}