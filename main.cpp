

#include <opencv2/core.hpp>

#include <iostream>

struct Vertex {
  double x_ = 0;
  double y_ = 0;
  double z_ = 0;

  double n_x_ = 0;
  double n_y_ = 0;
  double n_z_ = 0;
};

struct Ray {
  // Origin
  double x_ = 0;
  double y_ = 0;
  double z_ = 0;

  // Direction
  double d_x_ = 0;
  double d_y_ = 0;
  double d_z_ = 0;

  // Color
  double r_ = 0;
  double g_ = 0;
  double b_ = 0;
};

struct Triangle {
  Vertex point0_;
  Vertex point1_;
  Vertex point2_;
};

class Object {
public:
  Object() = default; // FIXME: Delete?
  virtual ~Object() = default;

  virtual double getVolume() = 0;

protected:
  std::vector<Triangle> triangles_;
};

class Camera {
public:
  Camera(const Vertex directed_center, const double distance_to_plane,
         const int width, const int height) : directed_center_(directed_center),
                                              distance_to_plane_(
                                                  distance_to_plane),
                                              width_(width), height_(height) {
  }

  void render();

//  void fireRay(const Ray ray,  );

private:
  Vertex directed_center_;
  double distance_to_plane_ = 0;

  const int width_;
  const int height_;
};

//void Camera::render() {
//  float angleUp = 0, angleRight = 0; // углы поворота от если бы мы смотрели по направлению 0.0.1 в радианах
//                                      //походу я изобрел сферические координаты и это альфа и тета
//  cv::Mat result_mat(height_, width_, CV_64FC3, cv::Scalar(0.0));
//
//  // Поиск координаты пикселя на плоскости в СК относительно центра камеры
//
//  for (size_t j = 0; j < height_; j++) {
//    for (size_t i = 0; i < width_; i++) {
//      float x = (2 * (i + 0.5) / width_ - 1) * width_ / (float) height_;
//      float y = -(2 * (j + 0.5) / height_ - 1);                           //y in (-1..1)
//      cv::Vec3f dir = cv::Vec3f(x, y, distance_to_plane_);            //стреляем луч из 0.0.0 по направлению
//
//      dir = cv::Vec3f(dir[1] * cos(angleRight) - dir[3] * sin(angleRight),
//                  dir[2] * cos(angleUp) - dir[3] * sin(angleUp),
//                  dir[1] * sin(angleRight) + dir[2] * sin(angleUp) +
//                  dir[3] * cos(angleUp)
//      );                                                              //поворачиваем его по желанию
//      cv::norm(dir);
//
//      ////////////////////////////////////////////////////////////////////////////////////
//      result_mat.at<double>(i, j) = cast_ray(cv::Vec3f(directed_center_.x_,directed_center_.y_,directed_center_.z_),
//                                             dir);
//    }
//  }
//}

class Scene {

};

int main() {

  // View
//  cv::namedWindow("Display window", cv::WINDOW_AUTOSIZE);
//  imshow("Display window", image);
//  cv::waitKey(0);
//  cv::imwrite("img.png", image);
  return 0;
}
