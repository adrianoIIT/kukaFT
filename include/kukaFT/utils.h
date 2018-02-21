/**
 \author Edoardo Lamon
        \file utils.h
    \brief Utilities library for planning and control
*/

#ifndef UTILS_H
#define UTILS_H

#include <cmath>
#include <vector>
#include <string>
#include <eigen3/Eigen/Dense>
#define LBR_MNJ 7		/*!< Number of data for a Cartesian vector */

#ifndef M_PI
#define M_PI 3.14159265358979323
#endif

using namespace std;

class utils
{
    public:

    static Eigen::MatrixXd calcKukaJacob(double jointPos[])
    {
        double C1,S1,S2,C2,S3,C3,S4,C4,S5,C5,S6,C6,S7,C7;
        double l1=0.310, l2=0.4, l3=0.39, l4=0.075;
        Eigen::MatrixXd J_kuka(FRI_CART_VEC,LBR_MNJ);

        C1=cos(jointPos[0]);
        S1=sin(jointPos[0]);
        C2=cos(jointPos[1]+M_PI/2);
        S2=sin(jointPos[1]+M_PI/2);
        C3=cos(jointPos[2]);
        S3=sin(jointPos[2]);
        C4=cos(jointPos[3]);
        S4=sin(jointPos[3]);
        C5=cos(jointPos[4]);
        S5=sin(jointPos[4]);
        C6=cos(jointPos[5]);
        S6=sin(jointPos[5]);
        C7=cos(jointPos[6]);
        S7=sin(jointPos[6]);

    //    J_kuka.row(0) << - l3*(S4*(C1*S3 + C3*S1*S2) + C2*C4*S1) - l4*(C6*(S4*(C1*S3 + C3*S1*S2) + C2*C4*S1) - S6*(C5*(C4*(C1*S3 + C3*S1*S2) - C2*S1*S4) + S5*(C1*C3 - S1*S2*S3))) - l2*C2*S1, -C1*(l3*(C4*S2 - C2*C3*S4) + l2*S2 + l4*(S6*(C5*(S2*S4 + C2*C3*C4) - C2*S3*S5) + C6*(C4*S2 - C2*C3*S4))) , C2*S1*(l3*(C4*S2 - C2*C3*S4) + l2*S2 + l4*(S6*(C5*(S2*S4 + C2*C3*C4) - C2*S3*S5) + C6*(C4*S2 - C2*C3*S4))) - S2*(l3*(S4*(C1*S3 + C3*S1*S2) + C2*C4*S1) + l4*(C6*(S4*(C1*S3 + C3*S1*S2) + C2*C4*S1) - S6*(C5*(C4*(C1*S3 + C3*S1*S2) - C2*S1*S4) + S5*(C1*C3 - S1*S2*S3))) + l2*C2*S1) , (C1*C3 - S1*S2*S3)*(l3*(C4*S2 - C2*C3*S4) + l4*(S6*(C5*(S2*S4 + C2*C3*C4) - C2*S3*S5) + C6*(C4*S2 - C2*C3*S4))) - C2*S3*(l3*(S4*(C1*S3 + C3*S1*S2) + C2*C4*S1) + l4*(C6*(S4*(C1*S3 + C3*S1*S2) + C2*C4*S1) - S6*(C5*(C4*(C1*S3 + C3*S1*S2) - C2*S1*S4) + S5*(C1*C3 - S1*S2*S3)))) , (S4*(C1*S3 + C3*S1*S2) + C2*C4*S1)*(l3*(C4*S2 - C2*C3*S4) + l4*(S6*(C5*(S2*S4 + C2*C3*C4) - C2*S3*S5) + C6*(C4*S2 - C2*C3*S4))) - (l3*(S4*(C1*S3 + C3*S1*S2) + C2*C4*S1) + l4*(C6*(S4*(C1*S3 + C3*S1*S2) + C2*C4*S1) - S6*(C5*(C4*(C1*S3 + C3*S1*S2) - C2*S1*S4) + S5*(C1*C3 - S1*S2*S3))))*(C4*S2 - C2*C3*S4) , l4*(S6*(C5*(S2*S4 + C2*C3*C4) - C2*S3*S5) + C6*(C4*S2 - C2*C3*S4))*(S5*(C4*(C1*S3 + C3*S1*S2) - C2*S1*S4) - C5*(C1*C3 - S1*S2*S3)) + l4*(C6*(S4*(C1*S3 + C3*S1*S2) + C2*C4*S1) - S6*(C5*(C4*(C1*S3 + C3*S1*S2) - C2*S1*S4) + S5*(C1*C3 - S1*S2*S3)))*(S5*(S2*S4 + C2*C3*C4) + C2*C5*S3) , 0;
    //    J_kuka.row(1) << l2*C1*C2 - l4*(C6*(S4*(S1*S3 - C1*C3*S2) - C1*C2*C4) - S6*(C5*(C4*(S1*S3 - C1*C3*S2) + C1*C2*S4) + S5*(C3*S1 + C1*S2*S3))) - l3*(S4*(S1*S3 - C1*C3*S2) - C1*C2*C4), -S1*(l3*(C4*S2 - C2*C3*S4) + l2*S2 + l4*(S6*(C5*(S2*S4 + C2*C3*C4) - C2*S3*S5) + C6*(C4*S2 - C2*C3*S4))), - S2*(l3*(S4*(S1*S3 - C1*C3*S2) - C1*C2*C4) + l4*(C6*(S4*(S1*S3 - C1*C3*S2) - C1*C2*C4) - S6*(C5*(C4*(S1*S3 - C1*C3*S2) + C1*C2*S4) + S5*(C3*S1 + C1*S2*S3))) - l2*C1*C2) - C1*C2*(l3*(C4*S2 - C2*C3*S4) + l2*S2 + l4*(S6*(C5*(S2*S4 + C2*C3*C4) - C2*S3*S5) + C6*(C4*S2 - C2*C3*S4))), (C3*S1 + C1*S2*S3)*(l3*(C4*S2 - C2*C3*S4) + l4*(S6*(C5*(S2*S4 + C2*C3*C4) - C2*S3*S5) + C6*(C4*S2 - C2*C3*S4))) - C2*S3*(l3*(S4*(S1*S3 - C1*C3*S2) - C1*C2*C4) + l4*(C6*(S4*(S1*S3 - C1*C3*S2) - C1*C2*C4) - S6*(C5*(C4*(S1*S3 - C1*C3*S2) + C1*C2*S4) + S5*(C3*S1 + C1*S2*S3)))), (l3*(C4*S2 - C2*C3*S4) + l4*(S6*(C5*(S2*S4 + C2*C3*C4) - C2*S3*S5) + C6*(C4*S2 - C2*C3*S4)))*(S4*(S1*S3 - C1*C3*S2) - C1*C2*C4) - (C4*S2 - C2*C3*S4)*(l3*(S4*(S1*S3 - C1*C3*S2) - C1*C2*C4) + l4*(C6*(S4*(S1*S3 - C1*C3*S2) - C1*C2*C4) - S6*(C5*(C4*(S1*S3 - C1*C3*S2) + C1*C2*S4) + S5*(C3*S1 + C1*S2*S3)))), l4*(S6*(C5*(S2*S4 + C2*C3*C4) - C2*S3*S5) + C6*(C4*S2 - C2*C3*S4))*(S5*(C4*(S1*S3 - C1*C3*S2) + C1*C2*S4) - C5*(C3*S1 + C1*S2*S3)) + l4*(C6*(S4*(S1*S3 - C1*C3*S2) - C1*C2*C4) - S6*(C5*(C4*(S1*S3 - C1*C3*S2) + C1*C2*S4) + S5*(C3*S1 + C1*S2*S3)))*(S5*(S2*S4 + C2*C3*C4) + C2*C5*S3), 0;
    //    J_kuka.row(2) << 0 , S1*(l3*(S4*(C1*S3 + C3*S1*S2) + C2*C4*S1) + l4*(C6*(S4*(C1*S3 + C3*S1*S2) + C2*C4*S1) - S6*(C5*(C4*(C1*S3 + C3*S1*S2) - C2*S1*S4) + S5*(C1*C3 - S1*S2*S3))) + l2*C2*S1) - C1*(l3*(S4*(S1*S3 - C1*C3*S2) - C1*C2*C4) + l4*(C6*(S4*(S1*S3 - C1*C3*S2) - C1*C2*C4) - S6*(C5*(C4*(S1*S3 - C1*C3*S2) + C1*C2*S4) + S5*(C3*S1 + C1*S2*S3))) - l2*C1*C2), C2*S1*(l3*(S4*(S1*S3 - C1*C3*S2) - C1*C2*C4) + l4*(C6*(S4*(S1*S3 - C1*C3*S2) - C1*C2*C4) - S6*(C5*(C4*(S1*S3 - C1*C3*S2) + C1*C2*S4) + S5*(C3*S1 + C1*S2*S3))) - l2*C1*C2) + C1*C2*(l3*(S4*(C1*S3 + C3*S1*S2) + C2*C4*S1) + l4*(C6*(S4*(C1*S3 + C3*S1*S2) + C2*C4*S1) - S6*(C5*(C4*(C1*S3 + C3*S1*S2) - C2*S1*S4) + S5*(C1*C3 - S1*S2*S3))) + l2*C2*S1), (C1*C3 - S1*S2*S3)*(l3*(S4*(S1*S3 - C1*C3*S2) - C1*C2*C4) + l4*(C6*(S4*(S1*S3 - C1*C3*S2) - C1*C2*C4) - S6*(C5*(C4*(S1*S3 - C1*C3*S2) + C1*C2*S4) + S5*(C3*S1 + C1*S2*S3)))) - (l3*(S4*(C1*S3 + C3*S1*S2) + C2*C4*S1) + l4*(C6*(S4*(C1*S3 + C3*S1*S2) + C2*C4*S1) - S6*(C5*(C4*(C1*S3 + C3*S1*S2) - C2*S1*S4) + S5*(C1*C3 - S1*S2*S3))))*(C3*S1 + C1*S2*S3), (S4*(C1*S3 + C3*S1*S2) + C2*C4*S1)*(l3*(S4*(S1*S3 - C1*C3*S2) - C1*C2*C4) + l4*(C6*(S4*(S1*S3 - C1*C3*S2) - C1*C2*C4) - S6*(C5*(C4*(S1*S3 - C1*C3*S2) + C1*C2*S4) + S5*(C3*S1 + C1*S2*S3)))) - (l3*(S4*(C1*S3 + C3*S1*S2) + C2*C4*S1) + l4*(C6*(S4*(C1*S3 + C3*S1*S2) + C2*C4*S1) - S6*(C5*(C4*(C1*S3 + C3*S1*S2) - C2*S1*S4) + S5*(C1*C3 - S1*S2*S3))))*(S4*(S1*S3 - C1*C3*S2) - C1*C2*C4), l4*(C6*(S4*(S1*S3 - C1*C3*S2) - C1*C2*C4) - S6*(C5*(C4*(S1*S3 - C1*C3*S2) + C1*C2*S4) + S5*(C3*S1 + C1*S2*S3)))*(S5*(C4*(C1*S3 + C3*S1*S2) - C2*S1*S4) - C5*(C1*C3 - S1*S2*S3)) - l4*(C6*(S4*(C1*S3 + C3*S1*S2) + C2*C4*S1) - S6*(C5*(C4*(C1*S3 + C3*S1*S2) - C2*S1*S4) + S5*(C1*C3 - S1*S2*S3)))*(S5*(C4*(S1*S3 - C1*C3*S2) + C1*C2*S4) - C5*(C3*S1 + C1*S2*S3)), 0;
    //    J_kuka.row(3) << 0 , S1 , C1*C2 , - C3*S1 - C1*S2*S3, C1*C2*C4 - S4*(S1*S3 - C1*C3*S2), C5*(C3*S1 + C1*S2*S3) - S5*(C4*(S1*S3 - C1*C3*S2) + C1*C2*S4), S6*(C5*(C4*(S1*S3 - C1*C3*S2) + C1*C2*S4) + S5*(C3*S1 + C1*S2*S3)) - C6*(S4*(S1*S3 - C1*C3*S2) - C1*C2*C4);
    //    J_kuka.row(4) << 0 , -C1 , C2*S1 , C1*C3 - S1*S2*S3, S4*(C1*S3 + C3*S1*S2) + C2*C4*S1, S5*(C4*(C1*S3 + C3*S1*S2) - C2*S1*S4) - C5*(C1*C3 - S1*S2*S3), C6*(S4*(C1*S3 + C3*S1*S2) + C2*C4*S1) - S6*(C5*(C4*(C1*S3 + C3*S1*S2) - C2*S1*S4) + S5*(C1*C3 - S1*S2*S3));
    //    J_kuka.row(5) << 1 , 0 , S2 , C2*S3 , C4*S2-C2*C3*S4 , - S5*(S2*S4+C2*C3*C4) - C2*C5*S3, S6*(C5*(S2*S4 + C2*C3*C4) - C2*S3*S5) + C6*(C4*S2 - C2*C3*S4);
        J_kuka <<- l3*(S4*(C1*S3 + C3*S1*S2) + C2*C4*S1) - l4*(C6*(S4*(C1*S3 + C3*S1*S2) + C2*C4*S1) - S6*(C5*(C4*(C1*S3 + C3*S1*S2) - C2*S1*S4) + S5*(C1*C3 - S1*S2*S3))) - l2*C2*S1, -C1*(l3*(C4*S2 - C2*C3*S4) + l2*S2 + l4*(S6*(C5*(S2*S4 + C2*C3*C4) - C2*S3*S5) + C6*(C4*S2 - C2*C3*S4))), C2*S1*(l3*(C4*S2 - C2*C3*S4) + l2*S2 + l4*(S6*(C5*(S2*S4 + C2*C3*C4) - C2*S3*S5) + C6*(C4*S2 - C2*C3*S4))) - S2*(l3*(S4*(C1*S3 + C3*S1*S2) + C2*C4*S1) + l4*(C6*(S4*(C1*S3 + C3*S1*S2) + C2*C4*S1) - S6*(C5*(C4*(C1*S3 + C3*S1*S2) - C2*S1*S4) + S5*(C1*C3 - S1*S2*S3))) + l2*C2*S1), (C1*C3 - S1*S2*S3)*(l3*(C4*S2 - C2*C3*S4) + l4*(S6*(C5*(S2*S4 + C2*C3*C4) - C2*S3*S5) + C6*(C4*S2 - C2*C3*S4))) - C2*S3*(l3*(S4*(C1*S3 + C3*S1*S2) + C2*C4*S1) + l4*(C6*(S4*(C1*S3 + C3*S1*S2) + C2*C4*S1) - S6*(C5*(C4*(C1*S3 + C3*S1*S2) - C2*S1*S4) + S5*(C1*C3 - S1*S2*S3)))), (S4*(C1*S3 + C3*S1*S2) + C2*C4*S1)*(l3*(C4*S2 - C2*C3*S4) + l4*(S6*(C5*(S2*S4 + C2*C3*C4) - C2*S3*S5) + C6*(C4*S2 - C2*C3*S4))) - (l3*(S4*(C1*S3 + C3*S1*S2) + C2*C4*S1) + l4*(C6*(S4*(C1*S3 + C3*S1*S2) + C2*C4*S1) - S6*(C5*(C4*(C1*S3 + C3*S1*S2) - C2*S1*S4) + S5*(C1*C3 - S1*S2*S3))))*(C4*S2 - C2*C3*S4), l4*(S6*(C5*(S2*S4 + C2*C3*C4) - C2*S3*S5) + C6*(C4*S2 - C2*C3*S4))*(S5*(C4*(C1*S3 + C3*S1*S2) - C2*S1*S4) - C5*(C1*C3 - S1*S2*S3)) + l4*(C6*(S4*(C1*S3 + C3*S1*S2) + C2*C4*S1) - S6*(C5*(C4*(C1*S3 + C3*S1*S2) - C2*S1*S4) + S5*(C1*C3 - S1*S2*S3)))*(S5*(S2*S4 + C2*C3*C4) + C2*C5*S3), 0,
               l2*C1*C2 - l4*(C6*(S4*(S1*S3 - C1*C3*S2) - C1*C2*C4) - S6*(C5*(C4*(S1*S3 - C1*C3*S2) + C1*C2*S4) + S5*(C3*S1 + C1*S2*S3))) - l3*(S4*(S1*S3 - C1*C3*S2) - C1*C2*C4), -S1*(l3*(C4*S2 - C2*C3*S4) + l2*S2 + l4*(S6*(C5*(S2*S4 + C2*C3*C4) - C2*S3*S5) + C6*(C4*S2 - C2*C3*S4))), - S2*(l3*(S4*(S1*S3 - C1*C3*S2) - C1*C2*C4) + l4*(C6*(S4*(S1*S3 - C1*C3*S2) - C1*C2*C4) - S6*(C5*(C4*(S1*S3 - C1*C3*S2) + C1*C2*S4) + S5*(C3*S1 + C1*S2*S3))) - l2*C1*C2) - C1*C2*(l3*(C4*S2 - C2*C3*S4) + l2*S2 + l4*(S6*(C5*(S2*S4 + C2*C3*C4) - C2*S3*S5) + C6*(C4*S2 - C2*C3*S4))), (C3*S1 + C1*S2*S3)*(l3*(C4*S2 - C2*C3*S4) + l4*(S6*(C5*(S2*S4 + C2*C3*C4) - C2*S3*S5) + C6*(C4*S2 - C2*C3*S4))) - C2*S3*(l3*(S4*(S1*S3 - C1*C3*S2) - C1*C2*C4) + l4*(C6*(S4*(S1*S3 - C1*C3*S2) - C1*C2*C4) - S6*(C5*(C4*(S1*S3 - C1*C3*S2) + C1*C2*S4) + S5*(C3*S1 + C1*S2*S3)))), (l3*(C4*S2 - C2*C3*S4) + l4*(S6*(C5*(S2*S4 + C2*C3*C4) - C2*S3*S5) + C6*(C4*S2 - C2*C3*S4)))*(S4*(S1*S3 - C1*C3*S2) - C1*C2*C4) - (C4*S2 - C2*C3*S4)*(l3*(S4*(S1*S3 - C1*C3*S2) - C1*C2*C4) + l4*(C6*(S4*(S1*S3 - C1*C3*S2) - C1*C2*C4) - S6*(C5*(C4*(S1*S3 - C1*C3*S2) + C1*C2*S4) + S5*(C3*S1 + C1*S2*S3)))), l4*(S6*(C5*(S2*S4 + C2*C3*C4) - C2*S3*S5) + C6*(C4*S2 - C2*C3*S4))*(S5*(C4*(S1*S3 - C1*C3*S2) + C1*C2*S4) - C5*(C3*S1 + C1*S2*S3)) + l4*(C6*(S4*(S1*S3 - C1*C3*S2) - C1*C2*C4) - S6*(C5*(C4*(S1*S3 - C1*C3*S2) + C1*C2*S4) + S5*(C3*S1 + C1*S2*S3)))*(S5*(S2*S4 + C2*C3*C4) + C2*C5*S3), 0,
               0, S1*(l3*(S4*(C1*S3 + C3*S1*S2) + C2*C4*S1) + l4*(C6*(S4*(C1*S3 + C3*S1*S2) + C2*C4*S1) - S6*(C5*(C4*(C1*S3 + C3*S1*S2) - C2*S1*S4) + S5*(C1*C3 - S1*S2*S3))) + l2*C2*S1) - C1*(l3*(S4*(S1*S3 - C1*C3*S2) - C1*C2*C4) + l4*(C6*(S4*(S1*S3 - C1*C3*S2) - C1*C2*C4) - S6*(C5*(C4*(S1*S3 - C1*C3*S2) + C1*C2*S4) + S5*(C3*S1 + C1*S2*S3))) - l2*C1*C2), C2*S1*(l3*(S4*(S1*S3 - C1*C3*S2) - C1*C2*C4) + l4*(C6*(S4*(S1*S3 - C1*C3*S2) - C1*C2*C4) - S6*(C5*(C4*(S1*S3 - C1*C3*S2) + C1*C2*S4) + S5*(C3*S1 + C1*S2*S3))) - l2*C1*C2) + C1*C2*(l3*(S4*(C1*S3 + C3*S1*S2) + C2*C4*S1) + l4*(C6*(S4*(C1*S3 + C3*S1*S2) + C2*C4*S1) - S6*(C5*(C4*(C1*S3 + C3*S1*S2) - C2*S1*S4) + S5*(C1*C3 - S1*S2*S3))) + l2*C2*S1), (C1*C3 - S1*S2*S3)*(l3*(S4*(S1*S3 - C1*C3*S2) - C1*C2*C4) + l4*(C6*(S4*(S1*S3 - C1*C3*S2) - C1*C2*C4) - S6*(C5*(C4*(S1*S3 - C1*C3*S2) + C1*C2*S4) + S5*(C3*S1 + C1*S2*S3)))) - (l3*(S4*(C1*S3 + C3*S1*S2) + C2*C4*S1) + l4*(C6*(S4*(C1*S3 + C3*S1*S2) + C2*C4*S1) - S6*(C5*(C4*(C1*S3 + C3*S1*S2) - C2*S1*S4) + S5*(C1*C3 - S1*S2*S3))))*(C3*S1 + C1*S2*S3), (S4*(C1*S3 + C3*S1*S2) + C2*C4*S1)*(l3*(S4*(S1*S3 - C1*C3*S2) - C1*C2*C4) + l4*(C6*(S4*(S1*S3 - C1*C3*S2) - C1*C2*C4) - S6*(C5*(C4*(S1*S3 - C1*C3*S2) + C1*C2*S4) + S5*(C3*S1 + C1*S2*S3)))) - (l3*(S4*(C1*S3 + C3*S1*S2) + C2*C4*S1) + l4*(C6*(S4*(C1*S3 + C3*S1*S2) + C2*C4*S1) - S6*(C5*(C4*(C1*S3 + C3*S1*S2) - C2*S1*S4) + S5*(C1*C3 - S1*S2*S3))))*(S4*(S1*S3 - C1*C3*S2) - C1*C2*C4), l4*(C6*(S4*(S1*S3 - C1*C3*S2) - C1*C2*C4) - S6*(C5*(C4*(S1*S3 - C1*C3*S2) + C1*C2*S4) + S5*(C3*S1 + C1*S2*S3)))*(S5*(C4*(C1*S3 + C3*S1*S2) - C2*S1*S4) - C5*(C1*C3 - S1*S2*S3)) - l4*(C6*(S4*(C1*S3 + C3*S1*S2) + C2*C4*S1) - S6*(C5*(C4*(C1*S3 + C3*S1*S2) - C2*S1*S4) + S5*(C1*C3 - S1*S2*S3)))*(S5*(C4*(S1*S3 - C1*C3*S2) + C1*C2*S4) - C5*(C3*S1 + C1*S2*S3)), 0,
               0, S1, C1*C2, - C3*S1 - C1*S2*S3, C1*C2*C4 - S4*(S1*S3 - C1*C3*S2), C5*(C3*S1 + C1*S2*S3) - S5*(C4*(S1*S3 - C1*C3*S2) + C1*C2*S4), S6*(C5*(C4*(S1*S3 - C1*C3*S2) + C1*C2*S4) + S5*(C3*S1 + C1*S2*S3)) - C6*(S4*(S1*S3 - C1*C3*S2) - C1*C2*C4),
               0, -C1, C2*S1, C1*C3 - S1*S2*S3, S4*(C1*S3 + C3*S1*S2) + C2*C4*S1, S5*(C4*(C1*S3 + C3*S1*S2) - C2*S1*S4) - C5*(C1*C3 - S1*S2*S3), C6*(S4*(C1*S3 + C3*S1*S2) + C2*C4*S1) - S6*(C5*(C4*(C1*S3 + C3*S1*S2) - C2*S1*S4) + S5*(C1*C3 - S1*S2*S3)),
               1, 0, S2, C2*S3, C4*S2 - C2*C3*S4, - S5*(S2*S4 + C2*C3*C4) - C2*C5*S3, S6*(C5*(S2*S4 + C2*C3*C4) - C2*S3*S5) + C6*(C4*S2 - C2*C3*S4);

        return J_kuka;

    }    
        
    static Eigen::VectorXd getScrewErrorFromAffine(Eigen::Affine3d &sourceAffine, Eigen::Affine3d &targetAffine)
    {
        Eigen::VectorXd screw(FRI_CART_VEC);
        Eigen::Affine3d targetFromSource = sourceAffine.inverse() * targetAffine;
        Eigen::AngleAxisd angleAxis(targetFromSource.linear());
        Eigen::Vector3d w = angleAxis.angle()*angleAxis.axis();
        screw << targetFromSource.translation(), w;

        return screw;
    }
    
    static Eigen::VectorXd getScrewFromAffine(Eigen::Affine3d &targetAffine)
    {
        Eigen::VectorXd screw(FRI_CART_VEC);
        Eigen::AngleAxisd angleAxis(targetAffine.linear());
        Eigen::Vector3d w = angleAxis.angle()*angleAxis.axis();
        screw << targetAffine.translation(), w;

        return screw;
    }

    static Eigen::Affine3d getAffineFromScrew(Eigen::VectorXd screw)
    {
//        cout << "-----------------------START getAffineFromScrew-----------------------" << endl;
        float angle;
        Eigen::Vector3d axisAngle, axis;
        Eigen::Matrix3d rotation;
        Eigen::Affine3d screwAffine;
        Eigen::Vector3d nullVec;
        nullVec.setZero(3);
//        Eigen::VectorXf myscrew(FRI_CART_VEC);

//        myscrew = screw;
//        cout << "Input screw : " << screw.transpose() << endl;
        axisAngle = screw.segment(3,3);
//        cout << "axisAngle: \n" << axisAngle.transpose() << endl;
        angle = axisAngle.norm();
        if (axisAngle == nullVec)
        {
            axis = nullVec;
        }
        else
            axis = axisAngle.normalized();
//        cout << "Angle is " << angle << " and axis is \n" << axis << endl;
        rotation = Eigen::AngleAxisd(angle,axis);
//         if (abs(rotation(0,0))==1 || abs(rotation(1,1))==1 || abs(rotation(2,2))==1)
//         {
//            cout << "matrix cleaning ..." << endl;
//         }
        if (abs(rotation(0,0))==1)
        {
            rotation(0,1) = 0.;
            rotation(0,2) = 0.;
            rotation(1,0) = 0.;
            rotation(2,0) = 0.;
        }
        if (abs(rotation(1,1))==1)
        {
            rotation(0,1) = 0.;
            rotation(1,0) = 0.;
            rotation(1,2) = 0.;
            rotation(2,1) = 0.;
        }
        if (abs(rotation(2,2))==1)
        {
            rotation(2,1) = 0.;
            rotation(0,2) = 0.;
            rotation(1,2) = 0.;
            rotation(2,0) = 0.;
        }
        screwAffine.linear() = rotation;
        screwAffine.translation() = screw.head(3);

//        cout << "-----------------------END getAffineFromScrew-----------------------" << endl;
        return screwAffine;
    }

    static vector<float> computeCubicTrajectory(float &relativePos, float avgJntVel) /*float avgJntVel=0.174)*/
    {
        //float finalTime = abs(relativePos / avgJntVel);
        float finalTime = 5.;
    //    cout << "Planning time: " << finalTime << " seconds" << endl;

        vector<float> trajCoeffs(4);
        trajCoeffs[0] = 0; // initial position
        trajCoeffs[1] = 0; // initial velocity
        trajCoeffs[2] = (3*(relativePos) / (pow(finalTime, 2)));   /*(-3*(initialJntPos - finalJntPos) - (2*initialJntVel + finalJntVel)) / (pow(finalTime, 2));*/
        trajCoeffs[3] = (-2*(relativePos)) / (pow(finalTime, 3));  /*(2*(initialJntPos - finalJntPos) + (initialJntVel + finalJntVel)) / (pow(finalTime, 3));*/

        return trajCoeffs;
    }


    static vector<double> computeQuinticTrajectory(double &relativePos, float finalTime=5.)
    {
//         float finalTime = 5.;
    //    cout << "Planning time: " << finalTime << " seconds" << endl;

        vector<double> trajCoeffs(6);
	
        trajCoeffs[0] = 0; // initial position
        trajCoeffs[1] = 0; // initial velocity
        trajCoeffs[2] = 0; // intial acceleration/2
        trajCoeffs[3] = (10*(relativePos)) / (pow(finalTime, 3));
        trajCoeffs[4] = (-15*(relativePos)) / (pow(finalTime, 4));
        trajCoeffs[5] = (6*(relativePos)) / (pow(finalTime, 5));

        return trajCoeffs;
    }

    static double nextStepCubic(vector<double> &cubicCoeffs, float t)
    {
        return cubicCoeffs[3]*pow(t,3.) + cubicCoeffs[2]*pow(t,2.) + cubicCoeffs[1]*t + cubicCoeffs[0]; // static_cast<float>
    }

    static double nextStepQuintic(vector<double> &quinticCoeffs, float t)
    {
        return quinticCoeffs[5]*pow(t,5.) + quinticCoeffs[4]*pow(t,4.) + quinticCoeffs[3]*pow(t,3.) + quinticCoeffs[2]*pow(t,2.) + quinticCoeffs[1]*t + quinticCoeffs[0];
    }
    
    static Eigen::Affine3d nextStepCircular(Eigen::Affine3d &initAffine, double t, double period, double radius, Eigen::Vector3d &center)
    {
	double omega = (2*M_PI)/(period);
	double x = - radius*cos(omega*t) + center(0); //+ initAffine.translation()(0);
	double y = radius*sin(omega*t) + center(1); // + initAffine.translation()(1);
	
// 	Eigen::VectorXf nextScrew(6);
// 	nextScrew.setZero(6);
// 	nextScrew(0) = x;
// 	nextScrew(1) = y;
// 	nextScrew(2) = initAffine.translation()(2);
	
	Eigen::Affine3d nextAffine;
	nextAffine.linear() = initAffine.linear();
	nextAffine.translation()(0) = x;
	nextAffine.translation()(1) = y;
	nextAffine.translation()(2) = initAffine.translation()(2);
	
	return nextAffine;
    }
    
    static Eigen::VectorXd computeForceImpedenceControl(Eigen::VectorXd posError, Eigen::VectorXd velError, Eigen::MatrixXd stiffness, Eigen::MatrixXd damping)
    {
        return  stiffness * posError + damping * velError;
    }

};


#endif // UTILS_H

