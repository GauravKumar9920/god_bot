from launch import LaunchDescription
from launch.actions import (
    DeclareLaunchArgument,
    IncludeLaunchDescription,
    SetEnvironmentVariable,
    TimerAction,
)
from launch.substitutions import LaunchConfiguration
from launch_ros.actions import Node
from ament_index_python.packages import get_package_share_directory
from launch.launch_description_sources import PythonLaunchDescriptionSource
import os


def generate_launch_description():
    pkg = get_package_share_directory('h_robot_planner')

    # Args
    model_name = LaunchConfiguration('model_name')
    declare_model_name = DeclareLaunchArgument(
        'model_name', default_value='omni_bot',
        description='Name for the spawned model')

    # Always start the stock empty world (world name == "empty")
    gz_sim_launch = IncludeLaunchDescription(
        PythonLaunchDescriptionSource(
            os.path.join(get_package_share_directory('ros_gz_sim'), 'launch', 'gz_sim.launch.py')
        ),
        launch_arguments={
            'gz_args': ['empty.sdf -r']
        }.items()
    )

    # Ensure GZ can find our package models for any future model:// URIs
    set_res_path = SetEnvironmentVariable(
        name='GZ_SIM_RESOURCE_PATH',
        value=os.pathsep.join([
            os.path.join(pkg, 'models'),
            os.environ.get('GZ_SIM_RESOURCE_PATH', '')
        ])
    )

    # Spawn from FILE PATH *and* target the known world name explicitly.
    model_sdf = os.path.join(pkg, 'models', 'omni_bot', 'model.sdf')
    spawn = Node(
        package='ros_gz_sim', executable='create', name='spawn_omni_bot',
        arguments=[
            '-world', 'empty',            # target world explicitly
            '-file', model_sdf,           # file path is bulletproof
            '-name', model_name,
            '-x', '0', '-y', '0', '-z', '0.1'
        ],
        output='screen')

    # Delay spawn slightly so server is definitely up
    spawn_delayed = TimerAction(period=2.0, actions=[spawn])

    # Bridge /model/omni_bot/cmd_vel (ROS<->GZ)
    bridge_cmd_vel = Node(
        package='ros_gz_bridge', executable='parameter_bridge', name='bridge_cmd_vel',
        arguments=[
            '/model/omni_bot/cmd_vel@geometry_msgs/msg/Twist@gz.msgs.Twist',
        ],
        output='screen')

    bridge_delayed = TimerAction(period=2.5, actions=[bridge_cmd_vel])

    return LaunchDescription([
        declare_model_name,
        set_res_path,
        gz_sim_launch,
        spawn_delayed,
        bridge_delayed,
    ])