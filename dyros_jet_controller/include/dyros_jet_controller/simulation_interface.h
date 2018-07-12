#ifndef SIMULATION_INTERFACE_H
#define SIMULATION_INTERFACE_H

#include "control_base.h"
#include "math_type_define.h"
#include <std_msgs/Float32.h>
#include <sensor_msgs/JointState.h>
#include <geometry_msgs/WrenchStamped.h>
#include <sensor_msgs/Imu.h>

// Used data structures:

namespace dyros_jet_controller
{

class SimulationInterface : public ControlBase{
public:
  SimulationInterface(ros::NodeHandle &nh, double Hz); // constructor for initialize node
  virtual ~SimulationInterface() { vrepStop(); }

  virtual void update() override; // update controller based on readdevice
  virtual void compute() override; // compute algorithm and update all class object
  virtual void writeDevice() override; // publish to actuate devices
  virtual void wait() override;

private:  // CALLBACK
  void simulationTimeCallback(const std_msgs::Float32ConstPtr& msg);
  void simulationStepDoneCallback(const std_msgs::BoolConstPtr& msg);
  void jointCallback(const sensor_msgs::JointStateConstPtr& msg);
  void leftFTCallback(const geometry_msgs::WrenchStampedConstPtr& msg);
  void rightFTCallback(const geometry_msgs::WrenchStampedConstPtr& msg);
  void gyroCallback(const geometry_msgs::PointConstPtr& msg);
  void accelCallback(const geometry_msgs::PointConstPtr& msg);
  void comCallback(const geometry_msgs::PointConstPtr& msg);
  void rfootPosCallback(const geometry_msgs::PointConstPtr& msg);
  void lfootPosCallback(const geometry_msgs::PointConstPtr& msg);
  void basePosCallback(const geometry_msgs::PointConstPtr& msg);
  void rfootOriCallback(const geometry_msgs::PointConstPtr& msg);
  void lfootOriCallback(const geometry_msgs::PointConstPtr& msg);
  void baseOriCallback(const geometry_msgs::PointConstPtr& msg);


private:
  void vrepStart();
  void vrepStop();
  void vrepStepTrigger();
  void vrepEnableSyncMode();

private:

  ros::Publisher vrep_joint_set_pub_;
  ros::Publisher vrep_sim_start_pub_;
  ros::Publisher vrep_sim_stop_pub_;
  ros::Publisher vrep_sim_step_trigger_pub_;
  ros::Publisher vrep_sim_enable_syncmode_pub_;

  ros::Subscriber vrep_sim_step_done_sub_;

  sensor_msgs::JointState joint_set_msg_;


  bool simulation_running_;
  bool simulation_step_done_;
  float simulation_time_; // from v-rep simulation time

  ros::Rate rate_;

  ros::Subscriber vrep_sim_state_sub_;

  ros::Subscriber imu_sub_;
  ros::Subscriber joint_sub_;
  ros::Subscriber left_ft_sub_;
  ros::Subscriber right_ft_sub_;
  ros::Subscriber com_sub_;
  ros::Subscriber gyro_sub_;
  ros::Subscriber accel_sub_;
  ros::Subscriber rfoot_pos_;
  ros::Subscriber lfoot_pos_;
  ros::Subscriber base_pos_;
  ros::Subscriber rfoot_ori_;
  ros::Subscriber lfoot_ori_;
  ros::Subscriber base_ori_;

};

}

#endif
