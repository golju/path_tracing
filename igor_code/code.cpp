#define _USE_MATH_DEFINES
#define M_PI
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <map>
#include <fstream>

using namespace std;

struct vec3 {
    double x = 0;
    double y = 0;
    double z = 0;

    vec3() {}

    vec3(double new_x, double new_y, double new_z) {
        x = new_x;
        y = new_y;
        z = new_z;
    }

    vec3(const vec3& vector) {
        if (&vector != this) {
            this->x = vector.x;
            this->y = vector.y;
            this->z = vector.z;
        }
    }

    vec3 operator + (const vec3& vector) const {
        return vec3(x + vector.x, y + vector.y, z + vector.z);
    }

    vec3 operator - (const vec3& vector) const {
        return vec3(x - vector.x, y - vector.y, z - vector.z);
    }

    vec3 operator * (double d) const {
        return vec3(x * d, y * d, z * d);
    }

    vec3 operator / (double d) const {
        return vec3(x / d, y / d, z / d);
    }

    vec3 normalize() const {
        double len = sqrt(x*x + y*y + z*z);
        return vec3(x / len, y / len, z / len);
    }

    double length() const {
        return sqrt(x * x + y * y + z * z);
    }
};

// Type aliases for vec3
using point3 = vec3;   // 3D point
using color = vec3;    // RGB color

double dot(const vec3& a, const vec3& b) {
    return (a.x * b.x + a.y * b.y + a.z * b.z);
}

vec3 cross(const vec3& a, const vec3& b) {
    return vec3(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x);
}

struct Ray {
    vec3 origin;
    vec3 direction;
    map <int, double> bright_coefficient; //spec
    map <int, double> L;

    Ray() {}

    Ray(const vec3& origin, const vec3& direction) {
        this->origin = origin;
        this->direction = direction;
    }

    Ray(const vec3& new_origin, const vec3& new_direction, const map<int, double>& bright_coefs, const map<int, double>& L) {
        this->origin = origin;
        this->direction = direction;
        this->bright_coefficient = bright_coefs;
        this->L = L;
    }
};

struct Sphere {
    vec3 center;
    double radius = 0;
    int sphere_id = -1;

    Sphere() {}

    Sphere(const vec3& new_center, double new_radius) {
        center = new_center;
        radius = new_radius;
    }

    vec3 getNormal(const vec3& inter_point) const {
        return (inter_point - center) / radius;
    }

    bool intersect(const Ray& ray, double& t) const {
        const vec3 o = ray.origin;
        const vec3 d = ray.direction;
        const vec3 oc = o - center;
        const double b = 2 * dot(oc, d);
        const double c = dot(oc, oc) - radius * radius;
        double disc = b * b - 4 * c;
        if (disc < 1e-4) return false;
        disc = sqrt(disc);
        const double t0 = -b - disc;
        const double t1 = -b + disc;
        t = (t0 < t1) ? t0 : t1;
        return true;
    }
};

struct Material {
    //double kd = -1; // Diffuse reflection coef 1 = 100%
    map <int, double> spec_Kd_color; // int wavelength, Kd color

    Material() {}

    Material(const  map <int, double>& spec_Kd_color) {
        this->spec_Kd_color = spec_Kd_color;
    }
};

struct Triangle {
    vec3 v0, v1, v2; 
    int object_id = -1;
    Material material;

    Triangle() {}

    Triangle(const vec3& v0, const vec3& v1, const vec3& v2) {
        this->v0 = v0;
        this->v1 = v1;
        this->v2 = v2;
    }

    Triangle(const vec3& v0, const vec3& v1, const vec3& v2, const int& object_id) {
        this->v0 = v0;
        this->v1 = v1;
        this->v2 = v2;
        this->object_id = object_id;
    }

    Triangle(const vec3& v0, const vec3& v1, const vec3& v2, const int& object_id, const Material& material) {
        this->v0 = v0;
        this->v1 = v1;
        this->v2 = v2;
        this->object_id = object_id;
        this->material = material;
    }

    vec3 getNormal() const {
        vec3 normal = cross(v1-v0, v2-v0);
        normal = normal.normalize();
        return normal;
    }

    vec3 getNormalByObserver(const vec3& observer) const {
        vec3 normal = cross(v1 - v0, v2 - v0);
        normal = normal.normalize();

        if (dot(observer.normalize(), normal) < 0) {
            normal = cross(v2-v0, v1-v0);
            normal = normal.normalize();
        }

        return normal;
    }

    bool intersect(const Ray& ray, double& t) const {
        vec3 e1 = v1 - v0;
        vec3 e2 = v2 - v0;
        vec3 pvec = cross(ray.direction, e2);
        double det = dot(e1, pvec);

        if (det < 1e-8 && det > -1e-8) return 0;

        double inv_det = 1 / det;
        vec3 tvec = ray.origin - v0;
        double u = dot(tvec, pvec) * inv_det;
        if (u < 0 || u > 1) return 0;

        vec3 qvec = cross(tvec, e1);
        double v = dot(ray.direction, qvec) * inv_det;
        if (v < 0 || v + u > 1) return 0;

        t = dot(e2, qvec) * inv_det;
        return t > 1e-8;
    }
};

struct Light {
    vec3 position;
    double total_intensity = 0; // W/sr
    map <int, double> spec_intensity; // int wavelength, double intensity

    Light() {};

    Light(const vec3& new_position, double new_total_flux) {
        position = new_position;
        total_intensity = new_total_flux;
    }

    Light(const vec3& new_position, double new_total_flux, const  map <int, double>& new_spec_intensity) {
        position = new_position;
        total_intensity = new_total_flux;
        spec_intensity = new_spec_intensity;
    }
};

bool scene_intersect(const Ray& ray, const vector<Triangle>& triangles, vec3 &hit, vec3 &N, Material &material) {

    double triangle_dist = numeric_limits<double>::max();

    for (int i = 0; i < triangles.size(); i++) {
        double t;
        if (triangles[i].intersect(ray, t) && t < triangle_dist) {
            triangle_dist = t;
            hit = ray.origin + ray.direction * t;    
            N = triangles[i].getNormalByObserver(ray.origin - hit);
            material = triangles[i].material;
        }
    }

    return triangle_dist < numeric_limits<double>::max();
}

Ray cast_ray(Ray& ray, const vector<Triangle>& triangles, const vector<Light>& lights) {
    vec3 hit, N;
    Material material;

    if (!scene_intersect(ray, triangles, hit, N, material)) {
        return ray;
    } 

    for (auto& item : material.spec_Kd_color) {
        ray.bright_coefficient.find(item.first)->second = (ray.bright_coefficient.find(item.first)->second) * item.second;
    }

    for (int i = 0; i < lights.size(); i++) {
        vec3 light_dir = (lights[i].position - hit).normalize();
        double dist = (lights[i].position - hit).length();
        double cos_teta = dot(light_dir, N);
        if (cos_teta <= 0) {
            return ray;
        }

        for (auto& item : ray.L) {
            double E = ((lights[i].spec_intensity.find(item.first)->second) / (dist * dist)) * cos_teta;
            item.second = (E * ray.bright_coefficient.find(item.first)->second) / (double)M_PI;
        }
    }

    return ray;
}

Ray fill_wavelength(Ray ray, vector<Light> lights) {
    map <int, double> bright_coefs;
    map <int, double> L;

    for (int i = 0; i < lights.size(); i++) {
        for (auto& item : lights[i].spec_intensity){
            bright_coefs.insert(make_pair(item.first, 1));
            L.insert(make_pair(item.first, 0));
        }
    }

    ray.bright_coefficient = bright_coefs;
    ray.L = L;

    return ray;
}

bool check(map<int, double> L) {

    for (auto& item : L)
    {
        if (item.second != 0) {
            return true;
        }
    }
    
    return false;
}

void render(const vector<Triangle>& triangles, const vector<Light>& lights) {
    const int width = 1024; //Было 1024
    const int height = 768; //Было 768
    const double fov = (double)M_PI / 3.f;
    vector<Ray> framebuffer(width * height);

    #pragma omp parallel for
    for (long long i = 0; i < height; i++) {
        for (long long j = 0; j < width; j++) {
            double x = -(2 * (j + 0.5) / (double)width - 1) * tan(fov / 2.) * width / (double)height;
            double y = -(2 * (i + 0.5) / (double)height - 1) * tan(fov / 2.);
            vec3 direction = vec3(x, y, -1).normalize(); //Было -1

            Ray ray(vec3(250, 300, 500), direction); //100, 200, 100 //250, 300, 500
            ray = fill_wavelength(ray, lights);

            framebuffer[j + i * width] = cast_ray(ray, triangles, lights);
        }
    }

    ofstream fout("../data/results_etalon.txt");

    vector<int> wavelengths;
    wavelengths.push_back(400);
    wavelengths.push_back(500);
    wavelengths.push_back(600);
    wavelengths.push_back(700);

    for (int k = 0; k < wavelengths.size(); k++) {
        fout << "wavelength" << " " << wavelengths[k] << endl;
        for (long long i = 0; i < height; i++) {
            for (long long j = 0; j < width; j++) {
                if (j == 0) {
                    fout << framebuffer[j + i * width].L.find(wavelengths[k])->second;
                }
                else {
                    fout << " " << framebuffer[j + i * width].L.find(wavelengths[k])->second;
                }
            }
            fout << endl;
        }
        fout << "\n";
    }

    fout.close();
}

int main() {

    //cornel_box0.shp

    vector<vec3> points;
    vector<Triangle> triangles;
    vector<int> object_id;
    vector<Light> lights;

    //Хардкод, так как в файле нет спектральных свойств поверхностей
    map <int, double> spec_Kd_color_brs_0;
    spec_Kd_color_brs_0.insert(make_pair(400, 0.343));
    spec_Kd_color_brs_0.insert(make_pair(500, 0.747));
    spec_Kd_color_brs_0.insert(make_pair(600, 0.74));
    spec_Kd_color_brs_0.insert(make_pair(700, 0.737));
    Material brs_0 = Material(spec_Kd_color_brs_0);

    map <int, double> spec_Kd_color_brs_1;
    spec_Kd_color_brs_1.insert(make_pair(400, 0.092));
    spec_Kd_color_brs_1.insert(make_pair(500, 0.285));
    spec_Kd_color_brs_1.insert(make_pair(600, 0.16));
    spec_Kd_color_brs_1.insert(make_pair(700, 0.159));
    Material brs_1 = Material(spec_Kd_color_brs_1);

    map <int, double> spec_Kd_color_brs_2;
    spec_Kd_color_brs_2.insert(make_pair(400, 0.04));
    spec_Kd_color_brs_2.insert(make_pair(500, 0.058));
    spec_Kd_color_brs_2.insert(make_pair(600, 0.287));
    spec_Kd_color_brs_2.insert(make_pair(700, 0.642));
    Material brs_2 = Material(spec_Kd_color_brs_2);

    map <int, double> spec_Kd_color_brs_3;
    spec_Kd_color_brs_3.insert(make_pair(400, 0.343));
    spec_Kd_color_brs_3.insert(make_pair(500, 0.747));
    spec_Kd_color_brs_3.insert(make_pair(600, 0.74));
    spec_Kd_color_brs_3.insert(make_pair(700, 0.737));
    Material brs_3 = Material(spec_Kd_color_brs_3);

    map <int, double> spec_Kd_color_brs_4;
    spec_Kd_color_brs_4.insert(make_pair(400, 0.343));
    spec_Kd_color_brs_4.insert(make_pair(500, 0.747));
    spec_Kd_color_brs_4.insert(make_pair(600, 0.74));
    spec_Kd_color_brs_4.insert(make_pair(700, 0.737));
    Material brs_4 = Material(spec_Kd_color_brs_4);

    vector<Material> materials;
    materials.push_back(brs_0);
    materials.push_back(brs_1);
    materials.push_back(brs_2);
    materials.push_back(brs_3);
    materials.push_back(brs_4);

    string s;
    ifstream file("../data/cornel_box0.shp");
    int state = 0; //1 - define breps brs_ найден, 2 - Number of vertices найден, Number of triangles найден
    int points_size = 0;
    while (getline(file, s)) {

        if (s.find("Number of parts") != string::npos) {
            state = 0;
            continue;
        }

        if (s.find("Number of triangles") != string::npos) {
            state = 3;
            continue;
        }

        if (state == 3) {
            istringstream iss(s);
            vector<int> number_of_triangles;
            for (std::string s; iss >> s;)
                number_of_triangles.push_back(stoi(s));

            Triangle triangle = Triangle(   points[points_size + number_of_triangles[0]],
                                            points[points_size + number_of_triangles[1]],
                                            points[points_size + number_of_triangles[2]],
                                            object_id[object_id.size() - 1],
                                            materials[object_id[object_id.size() - 1]]
                                        );
            triangles.push_back(triangle);
            continue;
        }

        if (s.find("define breps brs_") != string::npos) {
            istringstream iss(s);
            string token;
            while (getline(iss, token, '_')) {}
            object_id.push_back(stoi(token));
            state = 1;
            continue;
        }

        if (s.find("Number of vertices") != string::npos) {
            state = 2;
            points_size = points.size();
            continue;
        }

        if (state == 2) {
            istringstream iss(s);
            vector<double> cords;
            for (std::string s; iss >> s; )
                cords.push_back(stod(s));
            vec3 new_point = vec3(cords[0], cords[1], cords[2]); // Было vec3(cords[0], cords[1], cords[2])
            points.push_back(new_point);
        }
    }

    file.close();

    cout << triangles.size() << endl;

    //Хардкод, так как в файле нет свойств и положения источника света
    int total_intensity = 1445872; // W/sr
    double Itotal = 2100;
    double Isim_400 = 0;
    double Isim_500 = 400;
    double Isim_600 = 780;
    double Isim_700 = 920;

    map <int, double> spec_intensity;
    
    spec_intensity.insert(make_pair(400, total_intensity * (Isim_400 / Itotal)));
    spec_intensity.insert(make_pair(500, total_intensity * (Isim_500 / Itotal)));
    spec_intensity.insert(make_pair(600, total_intensity * (Isim_600 / Itotal)));
    spec_intensity.insert(make_pair(700, total_intensity * (Isim_700 / Itotal)));

    lights.push_back(Light(vec3(278, 545, -279.5), 1445872, spec_intensity)); //Было Light(vec3(278, -279.5, 548.7), 1445872, spec_intensity)

    //278, 548.7, -279.5

    render(triangles, lights);

    return 0;

}