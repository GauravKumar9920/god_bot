cmake_minimum_required(VERSION 3.8)
project(h_robot_planner)

if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  add_compile_options(-Wall -Wextra -Wpedantic)
endif()

# find dependencies
find_package(ament_cmake REQUIRED)
find_package(rclcpp REQUIRED)
find_package(geometry_msgs REQUIRED)
find_package(nav_msgs REQUIRED)
find_package(tf2_ros REQUIRED)
find_package(ros_gz_sim REQUIRED)
find_package(ros_gz_bridge REQUIRED)
find_package(visualization_msgs REQUIRED)
# find_package(h_motion_planner REQUIRED)

install(DIRECTORY models DESTINATION share/${PROJECT_NAME})
install(DIRECTORY worlds DESTINATION share/${PROJECT_NAME})
install(DIRECTORY launch DESTINATION share/${PROJECT_NAME})

# Public headers (so other packages/nodes can include <h_robot_planner/...>)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

# Spline visualizer node (publishes nav_msgs/Path and visualization_msgs/Marker)
add_executable(spline_viz_node src/spline_viz_node.cpp)
ament_target_dependencies(spline_viz_node rclcpp geometry_msgs nav_msgs visualization_msgs)

target_include_directories(spline_viz_node PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

install(TARGETS spline_viz_node DESTINATION lib/${PROJECT_NAME})
install(DIRECTORY include/ DESTINATION include)

if(BUILD_TESTING)
  find_package(ament_lint_auto REQUIRED)
  # the following line skips the linter which checks for copyrights
  # comment the line when a copyright and license is added to all source files
  set(ament_cmake_copyright_FOUND TRUE)
  # the following line skips cpplint (only works in a git repo)
  # comment the line when this package is in a git repo and when
  # a copyright and license is added to all source files
  set(ament_cmake_cpplint_FOUND TRUE)
  ament_lint_auto_find_test_dependencies()
endif()

ament_package()
