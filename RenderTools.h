#ifndef RENDER_TOOLS_H
#define RENDER_TOOLS_H

#include <opencv2/core/matx.hpp>

#include <vector>
#include <cmath>
#include <map>

//==============================================================================
//================================ Free functions ==============================
//==============================================================================
cv::Vec3d get_normalized(const cv::Vec3d &vec);

double get_length(const cv::Vec3d &vec);

//==============================================================================
//================================ Point =======================================
//==============================================================================
struct Point {
  double x_ = 0;
  double y_ = 0;
  double z_ = 0;

  // FIXME: На будущее
//  double n_x_ = 0;
//  double n_y_ = 0;
//  double n_z_ = 0;
};

//==============================================================================
//================================ Ray =========================================
//==============================================================================
struct Ray {
  cv::Vec3d origin;
  cv::Vec3d direction;
  std::map<int, double> bright_coefs; //spec
  std::map<int, double> L;

  Ray() {
  }

  Ray(const cv::Vec3d &origin, const cv::Vec3d &direction) {
    this->origin = origin;
    this->direction = direction;
  }

  Ray(const cv::Vec3d &new_origin, const cv::Vec3d &new_direction,
      const std::map<int, double> &bright_coefs,
      const std::map<int, double> &L) {
    this->origin = origin;
    this->direction = direction;
    this->bright_coefs = bright_coefs;
    this->L = L;
  }
};

//==============================================================================
//================================ Material ====================================
//==============================================================================
// FIXME: Maybe smth is private
class Material {
public:
  //double kd = -1; // Diffuse reflection coef 1 = 100%
  std::map<int, double> spec_Kd_color; // int wavelength, Kd color

  Material() {
  }

  Material(const std::map<int, double> &spec_Kd_color) {
    this->spec_Kd_color = spec_Kd_color;
  }
};

//==============================================================================
//================================ Triangle ====================================
//==============================================================================
// FIXME: Maybe smth is private
class Triangle {
public:
  cv::Vec3d v0, v1, v2;
  int object_id = -1;
  Material material;

  Triangle() {
  }

  Triangle(const cv::Vec3d &v0, const cv::Vec3d &v1, const cv::Vec3d &v2) {
    this->v0 = v0;
    this->v1 = v1;
    this->v2 = v2;
  }

  Triangle(const cv::Vec3d &v0, const cv::Vec3d &v1, const cv::Vec3d &v2,
           const int &object_id) {
    this->v0 = v0;
    this->v1 = v1;
    this->v2 = v2;
    this->object_id = object_id;
  }

  Triangle(const cv::Vec3d &v0, const cv::Vec3d &v1, const cv::Vec3d &v2,
           const int &object_id,
           const Material &material) {
    this->v0 = v0;
    this->v1 = v1;
    this->v2 = v2;
    this->object_id = object_id;
    this->material = material;
  }

  cv::Vec3d getNormal() const {
    cv::Vec3d normal = (v1 - v0).cross(v2 - v0);
    normal = get_normalized(normal);
    return normal;
  }

  cv::Vec3d getNormalByObserver(const cv::Vec3d &observer) const {
    cv::Vec3d normal = (v1 - v0).cross(v2 - v0);

    normal = get_normalized(normal);

    if (get_normalized(observer).dot(normal) < 0) {
      normal = (v2 - v0).cross(v1 - v0);
      normal = get_normalized(normal);
    }

    return normal;
  }

  bool intersect(const Ray &ray, double &t) const {
    cv::Vec3d e1 = v1 - v0;
    cv::Vec3d e2 = v2 - v0;
    cv::Vec3d pvec = ray.direction.cross(e2);
    double det = e1.dot(pvec);

    if (det < 1e-8 && det > -1e-8)
      return 0;

    double inv_det = 1 / det;
    cv::Vec3d tvec = ray.origin - v0;
    double u = tvec.dot(pvec) * inv_det;
    if (u < 0 || u > 1)
      return 0;

    cv::Vec3d qvec = tvec.cross(e1);
    double v = ray.direction.dot(qvec) * inv_det;
    if (v < 0 || v + u > 1)
      return 0;

    t = e2.dot(qvec) * inv_det;
    return t > 1e-8;
  }
};

//==============================================================================
//================================ Light =======================================
//==============================================================================
struct Light {
  cv::Vec3d position;
  double total_intensity = 0; // W/sr
  std::map<int, double> spec_intensity; // int wavelength, double intensity

  Light() {
  };

  Light(const cv::Vec3d &new_position, double new_total_flux) {
    position = new_position;
    total_intensity = new_total_flux;
  }

  Light(const cv::Vec3d &new_position, double new_total_flux,
        const std::map<int, double> &new_spec_intensity) {
    position = new_position;
    total_intensity = new_total_flux;
    spec_intensity = new_spec_intensity;
  }
};

//==============================================================================
//================================ Camera ======================================
//==============================================================================
class Camera {
public:
  Camera();

  void getPicture();

  int getWidth() const;

  int getHeight() const;

  double getFov() const;

private:
  const int width = -1; //Было 1024
  const int height = -1; //Было 768
  const double fov = -1.0;

  const cv::Vec3d origin = {-1.0, -1.0, -1.0};
};

//==============================================================================
//================================ Scene =======================================
//==============================================================================
class Scene {
public:
  Scene();

  void render();

  Ray fireRay(Ray &ray);

  int load(const std::string &path_to_file);

private:
  bool
  intersect(const Ray &ray, cv::Vec3d &hit, cv::Vec3d &N, Material &material);

private:
  std::vector<cv::Vec3d> points;
  std::vector<Triangle> triangles;
  std::vector<int> object_id;
  std::vector<Light> lights;
  std::vector<Material> materials;
  std::vector<Camera> cameras;
};

#endif //RENDER_TOOLS_H
