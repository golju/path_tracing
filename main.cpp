#include "RenderTools.h"

#include <iostream>

int main() {
  std::string path_to_scene_description = "../data/cornel_box0.shp";

  Scene scene;
  if (int exit_code = scene.load(path_to_scene_description); exit_code != 0) {
    std::cout << "Scene: failed to load scene description: " << exit_code
              << std::endl;
  }

  scene.render();
  return 0;
}