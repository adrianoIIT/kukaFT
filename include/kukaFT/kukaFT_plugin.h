/*
 * Copyright (C) 2017 IIT-ADVR
 * Author:
 * email:
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>
*/


#ifndef kukaFT_PLUGIN_H_
#define kukaFT_PLUGIN_H_

#include <XCM/XBotControlPlugin.h>

//#include <eigen3/Eigen/Dense>



namespace XBotPlugin {

/**
 * @brief kukaFT XBot RT Plugin
 *
 **/
class kukaFT : public XBot::XBotControlPlugin
{

public:

    kukaFT();
    virtual bool init_control_plugin(XBot::Handle::Ptr handle);

    virtual bool close();

    virtual void on_start(double time);

    virtual void on_stop(double time);
    
    virtual ~kukaFT();

protected:

    virtual void control_loop(double time, double period);

private:

    XBot::RobotInterface::Ptr _robot;

    double _start_time;

    Eigen::VectorXd _q0, jntPos;

    XBot::MatLogger::Ptr _logger;
    
    std::string base_frame,end_effector;
    
    std::vector< std::vector<double> > quinticTrajCoeffs;
    
    Eigen::Vector6d screwError, desScrewError, screwErrorPrev, relativeScrewError,
                        desTwistError, twistError, relativeNextDesiredScrew, msrdCartScrew, desCartScrew;

    Eigen::Affine3d  msrdCartAffine, initialCartAffine, desCartAffine,
                        nextDesiredAffine, relativeNextDesiredAffine, errorAffine, relativeErrorAffine;
			
    Eigen::Matrix6d measToWorldRotation;

    double timeCounter, deltaT, smoothingFactor;
    
    bool firstLoopCmd;
    
    // Circular trajectory
    const double radius = 0.05;
    const double period = 5.;
    Eigen::Vector3d center;
};

}

#endif // kukaFT_PLUGIN_H_
