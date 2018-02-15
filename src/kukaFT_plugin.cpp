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

#include <kukaFT_plugin.h>

#define FRI_CART_VEC 6
#include <utils.h>

/* Specify that the class XBotPlugin::kukaFT is a XBot RT plugin with name "kukaFT" */
REGISTER_XBOT_PLUGIN_(XBotPlugin::kukaFT)

namespace XBotPlugin {

  

  
kukaFT::kukaFT(): quinticTrajCoeffs(6), jntPos(7){
	twistError.setZero(FRI_CART_VEC);
	screwError.setZero(FRI_CART_VEC);
	measToWorldRotation.setZero(FRI_CART_VEC, FRI_CART_VEC);
}

bool kukaFT::init_control_plugin(XBot::Handle::Ptr handle)
{
    /* This function is called outside the real time loop, so we can
     * allocate memory on the heap, print stuff, ...
     * The RT plugin will be executed only if this init function returns true. */


    /* Save robot to a private member. */
    _robot = handle->getRobotInterface();

    /* Initialize a logger which saves to the specified file. Remember that
     * the current date/time is always appended to the provided filename,
     * so that logs do not overwrite each other. */

    _logger = XBot::MatLogger::getLogger("/tmp/kukaFT_log");
    
    //quinticTrajCoeffs.reserve(FRI_CART_VEC);
    
//     Eigen::VectorXd screwError, desScrewError, screwErrorPrev, relativeScrewError,
//                         desTwistError, twistError, relativeNextDesiredScrew;

    
    return true;


    
}


void kukaFT::on_start(double time)
{
    /* This function is called on plugin start, i.e. when the start command
     * is sent over the plugin switch port (e.g. 'rosservice call /kukaFT_switch true').
     * Since this function is called within the real-time loop, you should not perform
     * operations that are not rt-safe. */

    /* Save the plugin starting time to a class member */
    _robot->getMotorPosition(_q0);
    Eigen::VectorXd nullStiffness(7);
    nullStiffness.setZero();
    _robot->setStiffness(nullStiffness);
    _robot->setDamping(nullStiffness);
    _robot->move();

    /* Save the robot starting config to a class member */
    _start_time = time;
     
    // Desired pose initialization
      //variable outside loops
       base_frame = _robot->model().chain("arm1").getBaseLinkName();
       end_effector = _robot->model().chain("arm1").getTipLinkName();
       deltaT = 0.001;
       _logger->add("deltaT",deltaT);
       timeCounter = 0;
       smoothingFactor = 0.1;
       firstLoopCmd = true;
      
		    _robot->model().getPose(end_effector,msrdCartAffine);
                    initialCartAffine = msrdCartAffine;
                    _logger->add("nextDesiredAffine", initialCartAffine.matrix());
		    Eigen::Matrix3d desRotation;
		    Eigen::Vector3d desTranslation;
		    desRotation.setIdentity();
// 		    desTranslation.setZero();
		    desTranslation << 0,-0.5,0;
                    desCartAffine.linear() =  desRotation;
                    desCartAffine.translation() = desTranslation;
                    desCartAffine = initialCartAffine * desCartAffine;
		    
		    
		    
                   Eigen::VectorXd relativeDesiredScrew = utils::getScrewErrorFromAffine(initialCartAffine, desCartAffine);
		   
                   for (int i = 0; i < FRI_CART_VEC; i++)
                   {
		     
                      quinticTrajCoeffs[i] = utils::computeQuinticTrajectory(relativeDesiredScrew(i));
		   }
		   
// 		   // Planning a circular trajectory (fix the center w.r.t. world frame)
//                    center(0) = initialCartAffine.translation()(0) + radius;
//                    center(1) = initialCartAffine.translation()(1);
//                    center(2) = initialCartAffine.translation()(2);
}



void kukaFT::on_stop(double time)
{
    /* This function is called on plugin stop, i.e. when the stop command
     * is sent over the plugin switch port (e.g. 'rosservice call /kukaFT_switch false').
     * Since this function is called within the real-time loop, you should not perform
     * operations that are not rt-safe. */
}


void kukaFT::control_loop(double time, double period)
{
    /* This function is called on every control loop from when the plugin is start until
     * it is stopped.
     * Since this function is called within the real-time loop, you should not perform
     * operations that are not rt-safe. */

    /* The following code checks if any command was received from the plugin standard port
     * (e.g. from ROS you can send commands with
     *         rosservice call /kukaFT_cmd "cmd: 'MY_COMMAND_1'"
     * If any command was received, the code inside the if statement is then executed. */
      //variable outside loops

      Eigen::VectorXd desiredTorque(7),desiredForce(6), calcTorque(7), gravity(7);
      Eigen::MatrixXd J,J_pinv;
      _robot->model().update();
      _robot->model().getPose(end_effector,msrdCartAffine);
      _logger->add("msrdCartAffine",msrdCartAffine.matrix());
      _robot->getJointPosition(jntPos);
      _robot->model().getJacobian(end_effector,J);
//       J = utils::calcKukaJacob(jntPos.data());
      
      msrdCartScrew = utils::getScrewFromAffine(msrdCartAffine);
      _logger->add("msrdCartScrew", msrdCartScrew);
      
      //Only for Gazebo simulation
      Eigen::VectorXd stiff_int(7);
      stiff_int << 0,0,0,0,0,0,0;
      _robot->setStiffness(stiff_int);
      _robot->setDamping(stiff_int);
      
       // Cartesian Impedence
        Eigen::MatrixXd myStiffness(FRI_CART_VEC,FRI_CART_VEC), myDamping(FRI_CART_VEC,FRI_CART_VEC);
        Eigen::VectorXd stiffnessVector(FRI_CART_VEC);
                        stiffnessVector << 800.,800.,800.,100.,100.,100.;
        myStiffness = stiffnessVector.asDiagonal();
        float dampingFactor = 0.7;
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigenSolver(myStiffness);
        Eigen::MatrixXd eigenvalDiag = eigenSolver.eigenvalues().asDiagonal();
        eigenvalDiag = eigenvalDiag.array().sqrt();
        myDamping = 2 * dampingFactor * eigenSolver.eigenvectors() * eigenvalDiag * eigenSolver.eigenvectors().transpose();
	
	timeCounter += deltaT;
	if (timeCounter <= 5.)
                        {
                             for (int i = 0; i < FRI_CART_VEC; i++)
                             {
                                 relativeNextDesiredScrew(i) = utils::nextStepQuintic(quinticTrajCoeffs[i], timeCounter);
                             }
                        }

		       /*desCartAffine = utils::nextStepCircular(initialCartAffine,timeCounter,period,radius,center);
                       _logger->add("desCartAffine", desCartAffine.matrix());
		       relativeNextDesiredScrew = utils::getScrewErrorFromAffine(initialCartAffine, desCartAffine);*/	  
		       
                       relativeNextDesiredAffine = utils::getAffineFromScrew(relativeNextDesiredScrew);
                       nextDesiredAffine = initialCartAffine * relativeNextDesiredAffine;
                       _logger->add("nextDesiredAffine", nextDesiredAffine.matrix());
		   
	               desCartScrew = utils::getScrewFromAffine(desCartAffine);
                       _logger->add("desCartScrew", desCartScrew);
                   

                        // Compute position and orientation error
                       relativeScrewError = utils::getScrewErrorFromAffine(msrdCartAffine, nextDesiredAffine);
                       _logger->add("relativeScrewError",relativeScrewError);  
                       measToWorldRotation.topLeftCorner(3,3) = msrdCartAffine.linear();
                       measToWorldRotation.bottomRightCorner(3,3) = msrdCartAffine.linear();
		       _logger->add("measToWorldRotation", measToWorldRotation);

                       desScrewError = measToWorldRotation * relativeScrewError;
                       _logger->add("desScrewError", desScrewError);
                       screwError = smoothingFactor * desScrewError + (1 - smoothingFactor)*screwError;
                       _logger->add("screwError", screwError);
		       
                       if (firstLoopCmd)                   
                       {
                           screwErrorPrev = desScrewError;
                       }
                       desTwistError = (desScrewError - screwErrorPrev) / deltaT;
                       _logger->add("desTwistError", desTwistError);
                       twistError = smoothingFactor * desTwistError + (1 - smoothingFactor)*twistError;
                       _logger->add("twistError", twistError);
                       screwErrorPrev = desScrewError;                  
                        
                        
      desiredForce = myStiffness * desScrewError + myDamping * twistError;
      _logger->add("desiredForce",desiredForce);
      
      _robot->model().computeGravityCompensation(gravity);
      _logger->add("gravity",gravity);
      
      calcTorque = J.transpose() * desiredForce;
      _logger->add("calcTorque",calcTorque);
      
      desiredTorque = J.transpose() * desiredForce + gravity;
      _logger->add("desiredTorque",desiredTorque);
      
      _robot->setEffortReference(desiredTorque);
      _robot->move();
      
}

bool kukaFT::close()
{
    /* This function is called exactly once, at the end of the experiment.
     * It can be used to do some clean-up, or to save logging data to disk. */

    /* Save logged data to disk */
    _logger->flush();

    return true;
}

kukaFT::~kukaFT()
{
  
}

}
