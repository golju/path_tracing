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
//================================ Ray =========================================
//==============================================================================
class Ray {
public:
  Ray() = default;
  Ray(const cv::Vec3d &origin, const cv::Vec3d &direction);
  Ray(const cv::Vec3d &new_origin, const cv::Vec3d &new_direction,
      std::map<std::string, double> bright_coefs,
      std::map<std::string, double> L);

public:
  cv::Vec3d origin;
  cv::Vec3d direction;
  std::map<std::string, double> bright_coefficient;
  std::map<std::string, double> L;
};

//==============================================================================
//================================ Material ====================================
//==============================================================================
class Material {
public:
  Material() = default;
  Material(std::map<std::string, double> rgb_Kd_color);

public:
  std::map<std::string, double> rgb_Kd_color;
};

//==============================================================================
//================================ Triangle ====================================
//==============================================================================
// FIXME: Maybe smth is private
class Triangle {
public:
  Triangle() = default;

  Triangle(const cv::Vec3d &v0, const cv::Vec3d &v1, const cv::Vec3d &v2);
  Triangle(const cv::Vec3d &v0, const cv::Vec3d &v1, const cv::Vec3d &v2,
           const int &object_id);
  Triangle(const cv::Vec3d &v0, const cv::Vec3d &v1, const cv::Vec3d &v2,
           const int &object_id,
           const Material *material);
  cv::Vec3d getNormalByObserver(const cv::Vec3d &observer) const;
  bool intersect(const Ray &ray, double &t) const;


public:
  const Material *material;
private:
  cv::Vec3d v0, v1, v2;
  int object_id = -1;
};

//==============================================================================
//================================ Light =======================================
//==============================================================================
struct Light {
public:
  Light() = default;
  Light(const cv::Vec3d &new_position, double new_total_flux);
  Light(const cv::Vec3d &new_position, double new_total_flux,
        std::map<std::string, double> new_spec_intensity);

public:
  cv::Vec3d position;
  double total_intensity = 0; // W/sr
  std::map<std::string, double> spec_intensity; // int wavelength, double intensity
};

//==============================================================================
//================================ Camera ======================================
//==============================================================================
class Camera {
public:
  Camera(int width, int height, double fov, const cv::Vec3d &origin);

  int getWidth() const;
  int getHeight() const;
  double getFov() const;

private:
  const int width = -1;
  const int height = -1;
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
  void setSavePath(const std::string &path_to_save_result);
  void setNewCamera(const Camera &camera);

private:
  bool
  intersect(const Ray &ray, cv::Vec3d &hit, cv::Vec3d &N, Material &material);

private:
  std::string path_to_save_result = "../data/results.txt";

  std::vector<cv::Vec3d> points;
  std::vector<Triangle> triangles;
  std::vector<int> object_id;
  std::vector<Light> lights;
  std::vector<Material> materials;
  std::vector<Camera> cameras;
};

#endif //RENDER_TOOLS_H