#include "dyros_jet_controller/dyros_jet_model.h"
#include "dyros_jet_controller/walking_controller.h"


namespace dyros_jet_controller
{

void WalkingController::linkMass()
{
  ///////////////////////////////////////////////////////////////////////////////////
  mass_total_ = 0;

  mass_body_(0) =   model_.getLinkMass(0);
  mass_body_(1) =   model_.getLinkMass(13);
  mass_body_(2) =   model_.getLinkMass(14);

  c_waist_[0] = model_.getLinkComPosition(0);
  c_waist_[1] = model_.getLinkComPosition(13);
  c_waist_[2] = model_.getLinkComPosition(14);
  for(unsigned int i = 1; i<7; i++)
  {
    mass_l_leg_(i-1) = model_.getLinkMass(i);
    c_l_leg_[i-1] = model_.getLinkComPosition(i);
  }
  for(unsigned int i = 7; i<13; i++)
  {
    mass_r_leg_(i-7) = model_.getLinkMass(i);
    c_r_leg_[i-7] = model_.getLinkComPosition(i);
  }
  for(unsigned int i = 15; i<22; i++)
  {
    mass_l_arm_(i-15) = model_.getLinkMass(i);
    c_l_arm_[i-15] = model_.getLinkComPosition(i);
  }
  for(unsigned int i = 22; i<29; i++)
  {
    mass_r_arm_(i-22) = model_.getLinkMass(i);
    c_r_arm_[i-22] = model_.getLinkComPosition(i);
  }

  for(unsigned int i=0; i<29; i++)
  {
    mass_total_ += model_.getLinkMass(i);
  }

  cout<<"mass_total: "<<mass_total_<<endl;

}

void WalkingController::linkInertia()
{
  ///////////////////////////////////////////////////////////////////////////////////
  inertia_total_.setZero();

  for(int i=0; i<29; i++)
  {
    Eigen::Isometry3d link_com_t_float;
    if(i == 0)
    {
      link_com_t_float.setIdentity();
    }
    else
    {
      link_com_t_float = model_.getCurrentLinkTransform(i-1);
    }

    link_com_t_float.translation() = DyrosMath::multiplyIsometry3dVector3d(link_com_t_float, model_.getLinkComPosition(i));

    inertia_link_float_[i] = inertiaTensorTransform(model_.getLinkInertia(i), model_.getLinkMass(i), link_com_t_float);
    inertia_total_ += inertia_link_float_[i];
  }
}

Eigen::Matrix3d WalkingController::inertiaTensorTransform(Eigen::Matrix3d local_inertia, double mass, Eigen::Isometry3d transformation)
{
  Eigen::Matrix3d inertia_prime;
  Eigen::Matrix3d R = transformation.linear();
  Eigen::Vector3d p = transformation.translation();

  double p_square = (p.transpose()*p);
  inertia_prime = R*local_inertia*R.transpose();
  inertia_prime += mass*((p_square)*Eigen::Matrix3d::Identity() - p*p.transpose());
  return inertia_prime;
}

void WalkingController::getComJacobian()
{

  Eigen::Matrix<double, 6, 6> j_leg_link_float[12];
  Eigen::Matrix<double, 6, 7> j_arm_link_float[14];
  Eigen::Matrix<double, 3, 6> j_leg_com_link_support[12];
  Eigen::Matrix<double, 3, 7> j_arm_com_link_support[14];
  Eigen::Isometry3d leg_link_transform[12];
  Eigen::Isometry3d arm_link_transform[14];
  Eigen::Matrix3d skew_c_leg;
  Eigen::Matrix3d skew_c_arm;

  Eigen::Matrix6d adjoint_leg_com[12];
  Eigen::Matrix6d adjoint_arm_com[12];


  for(int i=0; i<12; i++)
  {
    j_leg_link_float[i] = model_.getLegLinkJacobian(i);  //left first
    leg_link_transform[i] = model_.getCurrentLinkTransform(i);
    if(walking_tick_ <= 1)
    {
      cout<< "j_leg_link_float["<<i<<"]:\n"<< j_leg_link_float[i]<< endl;
      cout<< "leg_link_transform["<<i<<"]:\n"<< leg_link_transform[i].translation()<< endl;
    }
    
  }
  for(int i=0; i<14; i++)
  {
    j_arm_link_float[i] = model_.getArmLinkJacobian(i);  //left first
    arm_link_transform[i] = model_.getCurrentLinkTransform(i+14);
    if(walking_tick_ <= 1)
    {
      cout<< "j_arm_link_float["<<i<<"]:\n"<< j_arm_link_float[i]<< endl;
      cout<< "arm_link_transform["<<i<<"]:\n"<< arm_link_transform[i].translation()<< endl;
    }
  }

  j_rleg_com_total_support.setZero();  // in the body center coordinates
  j_lleg_com_total_support.setZero();
  j_rarm_com_total_support.setZero();
  j_larm_com_total_support.setZero();

  adjoint_support_.block<3, 3>(0, 0) = pelv_support_current_.linear();
  adjoint_support_.block<3, 3>(3, 3) = pelv_support_current_.linear();


  for(int i=0; i<12; i++)
  {
    if(i<6)
    {
      skew_c_leg = DyrosMath::skew(leg_link_transform[i].linear()*c_l_leg_[i]);
    }
    else
    {
      skew_c_leg = DyrosMath::skew(leg_link_transform[i].linear()*c_r_leg_[i-6]);
    }

    adjoint_leg_com[i].setIdentity();
    adjoint_leg_com[i].block<3, 3>(0, 3) = -skew_c_leg; //Spatial Translation Matrix from i_th link origin to i_th link com position with respect to base coordinate

    j_leg_com_link_support[i] = (adjoint_support_*adjoint_leg_com[i]*j_leg_link_float[i]).block<3, 6>(0, 0);

    if(i<6)
    {
      j_lleg_com_total_support += (mass_l_leg_(i)/mass_total_)*j_leg_com_link_support[i];
    }
    else
    {
      j_rleg_com_total_support += (mass_r_leg_(i-6)/mass_total_)*j_leg_com_link_support[i];
    }

    if(walking_tick_ <= 1)
    {
      cout<< "j_leg_com_link_support["<<i<<"]:\n"<< j_leg_com_link_support[i]<< endl;
    }
  }

  if(walking_tick_ <= 1)
  {
    cout<< "j_lleg_com_total_support_: \n"<< j_lleg_com_total_support<< endl;
  }

  if(walking_tick_ <= 1)
  {
    cout<< "j_rleg_com_total_support: \n"<< j_rleg_com_total_support<< endl;
  }


  for(int i=0; i<14; i++)
  {

    if(i<7)
    {
      skew_c_arm = DyrosMath::skew(arm_link_transform[i].linear()*c_l_arm_[i]);
    }
    else
    {
      skew_c_arm = DyrosMath::skew(arm_link_transform[i].linear()*c_r_arm_[i-7]);
    }
    adjoint_arm_com[i].setIdentity();
    adjoint_arm_com[i].block<3, 3>(0, 3) = -skew_c_arm; //Spatial Translation Matrix from i_th link origin to i_th link com position with respect to base coordinate


    j_arm_com_link_support[i] = (adjoint_support_*adjoint_arm_com[i]*j_arm_link_float[i]).block<3, 7>(0, 0);

    if(i<7)
    {
      j_larm_com_total_support += (mass_l_arm_(i)/mass_total_)*j_arm_com_link_support[i];
    }
    else
    {
      j_rarm_com_total_support += (mass_r_arm_(i-7)/mass_total_)*j_arm_com_link_support[i];
    }
  }



  double kc;
  double kp;
  double kd;
  double kf;
  double kw;
  double lambda;
  Eigen::Vector3d r_c1;
  Eigen::Matrix3d skew_r_c1;
  Eigen::Matrix3d skew_r2_r1;
  Eigen::Vector3d error_foot;
  Eigen::Vector4d error_foot_w;
  Eigen::Vector3d swing_foot_w;
  Eigen::Vector3d error_com;
  Eigen::Vector3d error_zmp;
  Eigen::Vector4d error_w;
  Eigen::Vector3d error_moment;

  double switch_l_ft;
  double switch_r_ft;
  kc = 180;     //sim desired: 150 sim real: 180 //com error gain
  kp = 30;      //sim desired: 30 sim real: 30  //zmp error gain
  kd = 0.000;     //dob gain
  kf = 180;    //sim desired: 150  sin real: 180  //foot position error gain 
  kw = 150;   //sim desired: 150, sim real: 150     //orientation error gain

  //kc = 300.0; kp = 45.0; kd = 0.005;  //gains for real robot
  //kf = 300.0; kw = 200.0;
  lambda = 0.000;
  error_zmp.setZero();

  if(estimator_flag_ == true)
  {
    error_com(0) = com_desired_(0) - X_hat_post_2_(0);
    error_com(1) = com_desired_(1) - X_hat_post_2_(1);
    error_com(2) = com_desired_(2) - com_support_current_(2);
  }
  else
  {
    error_com = com_desired_ - com_support_current_;
  }

  if(l_ft_(2)+r_ft_(2) < -250)
  {
    error_zmp(0) = zmp_desired_(0) - zmp_measured_(0);
    error_zmp(1) = zmp_desired_(1) - zmp_measured_(1);
    // if(walking_tick_ <= 10)
    // {
    //   cout <<"error_zmp: \n"<< error_zmp <<endl;
    // }
  }
  else
  {
    if(walking_tick_%100 == 0)
    {
      // cout<<"I'm flying"<<endl;
      // cout<<"l_ft_(2): "<<l_ft_(2)<<endl;
      // cout<<"r_ft_(2): "<<r_ft_(2)<<endl;
    }
  }


  if(l_ft_(2) > 10)
  {
    switch_l_ft = 1;
  }
  else
  {
    switch_l_ft = 0;
  }

  if(r_ft_(2) > 10)
  {
    switch_r_ft = 1;
  }
  else
  {
    switch_r_ft = 0;
  }






  disturbance_accel_old_ = disturbance_accel_;


  disturbance_accel_(0) = desired_u_dot_(0) - (switch_l_ft*(-l_ft_(4)) + switch_r_ft*(-r_ft_(4)))/(mass_total_*com_support_current_(2));
  disturbance_accel_(1) = desired_u_dot_(1) - (switch_l_ft*(l_ft_(3)) + switch_r_ft*(r_ft_(3)))/(mass_total_*com_support_current_(2));
  disturbance_accel_(2) = 0;

  //disturbance_accel_ = 0.3*disturbance_accel_ +0.7*disturbance_accel_old_;
  //cout<<"error_moment"<<error_moment<<endl;


  desired_u_old_ = desired_u_;
  desired_u_ = com_dot_desired_ + kc*(error_com) - kp*(error_zmp);
  //desired_u_ = com_dot_desired_ + kc*(error_com) - 1*(error_moment);
  //desired_u_ = com_dot_desired_ + kc*(error_com) - kd*(disturbance_accel_);

  desired_u_dot_ = (desired_u_ - desired_u_old_)*hz_;
  error_w = DyrosMath::rot2Axis(pelv_trajectory_support_.linear()*(pelv_support_current_.linear().transpose()));
  desired_w_ =  kw*(error_w.segment<3>(0)*error_w(3));
  if(walking_tick_ <= 1)
  {
    cout << "desired_u_: \n"<< desired_u_<<endl;
    cout << "desired_w_: \n"<< desired_w_<<endl;
    cout << "pelv_trajectory_support_.linear(): \n"<< pelv_trajectory_support_.linear()<<endl;
    cout << "pelv_support_current_.linear(): \n"<< pelv_support_current_.linear()<<endl;
  }
    
  if (foot_step_(current_step_num_, 6) == 1) //left support foot
  {


    r_c1 = com_support_current_ - lfoot_support_current_.translation();
    skew_r_c1 = DyrosMath::skew(r_c1);
    j1_ = adjoint_support_*current_leg_jacobian_l_;
    j2_ = adjoint_support_*current_leg_jacobian_r_;

    j_v1_ = j1_.block<3, 6>(0, 0);
    j_w1_ = j1_.block<3, 6>(3, 0);


    skew_r2_r1 = DyrosMath::skew(pelv_support_current_.linear()*(lfoot_float_current_.translation() - rfoot_float_current_.translation()));
    //Skew(FOOT.L_T.translation()-FOOT.R_T.translation(), skew_r2_r1);

    adjoint_21_.setIdentity();
    adjoint_21_.block<3, 3>(0, 3) = skew_r2_r1; //Spatial Translation Matrix from i_th link origin to i_th link com position with respect to base coordinate

    error_foot = rfoot_trajectory_support_.translation() - rfoot_support_current_.translation();
    error_foot_w = DyrosMath::rot2Axis(rfoot_trajectory_support_.linear() * (rfoot_support_current_.linear().transpose()));

    swing_foot_w.setZero();
    swing_foot_w(2) = rfoot_trajectory_dot_support_(5);

    // x2_d_dot_.segment<3>(0) = rfoot_trajectory_dot_support_.segment<3>(0) + kf*(error_foot);
    x2_d_dot_.segment<3>(0) = kf*(error_foot);
    x2_d_dot_.segment<3>(3) = swing_foot_w + kw*(error_foot_w.segment<3>(0)*error_foot_w(3));


    j_com_psem_ = -j_v1_ + skew_r_c1*j_w1_ + j_lleg_com_total_support + j_rleg_com_total_support*j2_.inverse()*adjoint_21_*j1_;

    desired_c_dot_psem_ = desired_u_ - j_rleg_com_total_support*j2_.inverse()*x2_d_dot_;



  }
  else //right support foot
  {
    r_c1 = com_support_current_ - rfoot_support_current_.translation();
    skew_r_c1 = DyrosMath::skew(r_c1);
    j1_ = adjoint_support_*current_leg_jacobian_r_;
    j2_ = adjoint_support_*current_leg_jacobian_l_;

    j_v1_ = j1_.block<3, 6>(0, 0);
    j_w1_ = j1_.block<3, 6>(3, 0);


    skew_r2_r1 = DyrosMath::skew(pelv_support_current_.linear()*(rfoot_float_current_.translation() - lfoot_float_current_.translation()));
    //Skew(FOOT.L_T.translation()-FOOT.R_T.translation(), skew_r2_r1);

    adjoint_21_.setIdentity();
    adjoint_21_.block<3, 3>(0, 3) = skew_r2_r1; //Spatial Translation Matrix from i_th link origin to i_th link com position with respect to base coordinate

    error_foot = lfoot_trajectory_support_.translation() - lfoot_support_current_.translation();
    error_foot_w = DyrosMath::rot2Axis(lfoot_trajectory_support_.linear() * (lfoot_support_current_.linear()).transpose());

    swing_foot_w.setZero();
    swing_foot_w(2) = lfoot_trajectory_dot_support_(5);

    // x2_d_dot_.segment<3>(0) = lfoot_trajectory_dot_support_.segment<3>(0) + kf*(error_foot);
    x2_d_dot_.segment<3>(0) = kf*(error_foot);
    x2_d_dot_.segment<3>(3) = swing_foot_w + kw*(error_foot_w.segment<3>(0)*error_foot_w(3));


    j_com_psem_ = -j_v1_ + skew_r_c1*j_w1_ + j_rleg_com_total_support + j_lleg_com_total_support*j2_.inverse()*adjoint_21_*j1_;

    desired_c_dot_psem_ = desired_u_ - j_lleg_com_total_support*j2_.inverse()*x2_d_dot_;

  }

  if(walking_tick_ <= 1)
  {
    cout << "j1_: \n"<< j1_<<endl;
    cout << "j2_: \n"<< j2_<<endl;
    cout << "j_com_psem_: \n"<< j_com_psem_<<endl;
    cout << "desired_c_dot_psem_: \n"<< desired_c_dot_psem_<<endl;
    cout << "error_foot: \n"<< error_foot<<endl;
    cout << "error_foot_w: \n"<< error_foot_w<<endl;
    cout << "x2_d_dot_: \n"<< x2_d_dot_<<endl;
    cout << "rfoot_trajectory_support_.linear(): \n"<< rfoot_trajectory_support_.linear()<<endl;
    cout << "rfoot_support_current_.linear(): \n"<< rfoot_support_current_.linear()<<endl;
    cout << "rfoot_trajectory_support_.linear() * (rfoot_support_current_.linear().transpose()): \n"<< rfoot_trajectory_support_.linear() * (rfoot_support_current_.linear().transpose())<<endl;
  }


}



void WalkingController::computeComJacobianControl(Eigen::Vector12d &desired_leg_q_dot)
{
  double lamb = 0.001;

  if (foot_step_(current_step_num_, 6) == 1) // left support foot
  {
    j_total_.block<3, 6>(0, 0) = j_com_psem_;
    j_total_.block<3, 6>(3, 0) = -j_w1_;

    c_total_.segment<3>(0) = desired_c_dot_psem_;
    c_total_.segment<3>(3) = desired_w_;


    desired_leg_q_dot.segment<6>(0) = j_total_.inverse()*c_total_;  //left
    desired_leg_q_dot.segment<6>(6) = j2_.inverse()*(x2_d_dot_ + adjoint_21_*j1_*desired_leg_q_dot.segment<6>(0));
  }
  else //right support foot
  {

    j_total_.block<3, 6>(0, 0) = j_com_psem_;
    j_total_.block<3, 6>(3, 0) = -j_w1_;

    c_total_.segment<3>(0) = desired_c_dot_psem_;
    c_total_.segment<3>(3) = desired_w_;


    desired_leg_q_dot.segment<6>(6) = j_total_.inverse()*c_total_; //right
    desired_leg_q_dot.segment<6>(0) = j2_.inverse()*(x2_d_dot_ + adjoint_21_*j1_*desired_leg_q_dot.segment<6>(6));
  }


}
}
