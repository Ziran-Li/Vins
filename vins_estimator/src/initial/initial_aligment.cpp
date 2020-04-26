#include "initial_alignment.h"
//陀螺仪偏置标定
//先用相机计算出从bk到bk+1下的旋转，然后用ＩＭＵ计算从bk+1到bk下的旋转
void solveGyroscopeBias(map<double, ImageFrame> &all_image_frame, Vector3d* Bgs)//整个方法的计算模型由论文给出
{
    Matrix3d A;
    Vector3d b;
    Vector3d delta_bg;
    A.setZero();
    b.setZero();
    map<double, ImageFrame>::iterator frame_i;
    map<double, ImageFrame>::iterator frame_j;
    for (frame_i = all_image_frame.begin(); next(frame_i) != all_image_frame.end(); frame_i++)
    {
        frame_j = next(frame_i);
        MatrixXd tmp_A(3, 3);
        tmp_A.setZero();
        VectorXd tmp_b(3);
        tmp_b.setZero();
        Eigen::Quaterniond q_ij(frame_i->second.R.transpose() * frame_j->second.R);//这个地方感觉写反了，应该是i帧的旋转乘j帧旋转的逆
        tmp_A = frame_j->second.pre_integration->jacobian.template block<3, 3>(O_R, O_BG);
        tmp_b = 2 * (frame_j->second.pre_integration->delta_q.inverse() * q_ij).vec();//某个特征值对应的特征向量
        A += tmp_A.transpose() * tmp_A;//其实就是处理AT*Ax=AT*b的问题，一般最小二乘直接处理Ａx=b的问题，即Ａx-b=0即可
        b += tmp_A.transpose() * tmp_b;//这里使用ＬＤＬＴ方法,因为AT*A一定是可逆的，两边同时乘以其逆就可以

    }//要想使用ldlt的方法，首先要保证Ａ矩阵可逆
    delta_bg = A.ldlt().solve(b);//使用ＬＴＬＤ方法求解最小二乘问题
    ROS_WARN_STREAM("gyroscope bias initial calibration " << delta_bg.transpose());

    for (int i = 0; i <= WINDOW_SIZE; i++)
        Bgs[i] += delta_bg;//将陀螺仪偏置存到Bgs[],最后再重新计算一次预积分

    for (frame_i = all_image_frame.begin(); next(frame_i) != all_image_frame.end( ); frame_i++)
    {
        frame_j = next(frame_i);
        frame_j->second.pre_integration->repropagate(Vector3d::Zero(), Bgs[0]);
    }
}


MatrixXd TangentBasis(Vector3d &g0)//切线基础
{
    Vector3d b, c;
    Vector3d a = g0.normalized();//归一化
    Vector3d tmp(0, 0, 1);
    if(a == tmp)
        tmp << 1, 0, 0;
    b = (tmp - a * (a.transpose() * tmp)).normalized();
    c = a.cross(b);//叉乘，c的方向应该是重力要迭代的方向
    MatrixXd bc(3, 2);
    bc.block<3, 1>(0, 0) = b;
    bc.block<3, 1>(0, 1) = c;
    return bc;
}
//g0 = ag1 + w1 * b1 + w2 * b2 
void RefineGravity(map<double, ImageFrame> &all_image_frame, Vector3d &g, VectorXd &x)//对求出来的重力进行进一步优化
{
    Vector3d g0 = g.normalized() * G.norm();//g0应该是要返回的值，G.norm()应该是９．８，此处不断迭代的是g0的方向
    Vector3d lx, ly;
    //VectorXd x;
    int all_frame_count = all_image_frame.size();
    int n_state = all_frame_count * 3 + 2 + 1;//此处的每帧状态向量包括，速度，重力加速度，和尺度；速度是三维的，重力加速度是２维的，尺度是１维的
//此处重力加速度修正的时候，重力加速度是２维的，重力加速度计算的时候是３维的
    MatrixXd A{n_state, n_state};
    A.setZero();
    VectorXd b{n_state};
    b.setZero();

    map<double, ImageFrame>::iterator frame_i;
    map<double, ImageFrame>::iterator frame_j;
    for(int k = 0; k < 4; k++)
    {
        MatrixXd lxly(3, 2);
        lxly = TangentBasis(g0);//切线基础，返回的是３＊２的矩阵
        int i = 0;
        for (frame_i = all_image_frame.begin(); next(frame_i) != all_image_frame.end(); frame_i++, i++)
        {
            frame_j = next(frame_i);//i,j是相邻的两帧

            MatrixXd tmp_A(6, 9);//重力修正的时候Ｈ矩阵是６＊９的，计算的时候是６＊１０
            tmp_A.setZero();
            VectorXd tmp_b(6);//b矩阵无论是计算还是修正，都是６＊１
            tmp_b.setZero();

            double dt = frame_j->second.pre_integration->sum_dt;

　　　　　　//修正时候的b矩阵与Ｈ矩阵，与计算时候是不一样的，Ｈ矩阵多了b向量，这个b是重力的分量
            tmp_A.block<3, 3>(0, 0) = -dt * Matrix3d::Identity();
            tmp_A.block<3, 2>(0, 6) = frame_i->second.R.transpose() * dt * dt / 2 * Matrix3d::Identity() * lxly;//与之前计算重力时，有所不同
            tmp_A.block<3, 1>(0, 8) = frame_i->second.R.transpose() * (frame_j->second.T - frame_i->second.T) / 100.0;     
            tmp_b.block<3, 1>(0, 0) = frame_j->second.pre_integration->delta_p + frame_i->second.R.transpose() * frame_j->second.R * TIC[0] - TIC[0] - frame_i->second.R.transpose() * dt * dt / 2 * g0;
　　　　　　　　//与之前计算重力时，有所不同
            tmp_A.block<3, 3>(3, 0) = -Matrix3d::Identity();
            tmp_A.block<3, 3>(3, 3) = frame_i->second.R.transpose() * frame_j->second.R;
            tmp_A.block<3, 2>(3, 6) = frame_i->second.R.transpose() * dt * Matrix3d::Identity() * lxly;//与之前计算重力时，有所不同
            tmp_b.block<3, 1>(3, 0) = frame_j->second.pre_integration->delta_v - frame_i->second.R.transpose() * dt * Matrix3d::Identity() * g0;//与之前计算重力时，有所不同


            Matrix<double, 6, 6> cov_inv = Matrix<double, 6, 6>::Zero();
            //cov.block<6, 6>(0, 0) = IMU_cov[i + 1];
            //MatrixXd cov_inv = cov.inverse();
            cov_inv.setIdentity();

            MatrixXd r_A = tmp_A.transpose() * cov_inv * tmp_A;
            VectorXd r_b = tmp_A.transpose() * cov_inv * tmp_b;

            A.block<6, 6>(i * 3, i * 3) += r_A.topLeftCorner<6, 6>();
            b.segment<6>(i * 3) += r_b.head<6>();

            A.bottomRightCorner<3, 3>() += r_A.bottomRightCorner<3, 3>();
            b.tail<3>() += r_b.tail<3>();

            A.block<6, 3>(i * 3, n_state - 3) += r_A.topRightCorner<6, 3>();
            A.block<3, 6>(n_state - 3, i * 3) += r_A.bottomLeftCorner<3, 6>();
        }
            A = A * 1000.0;
            b = b * 1000.0;
            x = A.ldlt().solve(b);//求状态向量，这个状态向量包括速度、重力加速度和尺度；速度是两帧之间的速度
            VectorXd dg = x.segment<2>(n_state - 3);//我感觉是减６，不是减３
            g0 = (g0 + lxly * dg).normalized() * G.norm();
            //double s = x(n_state - 1);
    }   
    g = g0;
}
//LinearAlignment计算尺度，重力加速度和速度，论文中给出的是相邻两个速度的模型，映射到整个n+1个速度模型中，
//A矩阵一定是一个对称矩阵，代码中的Ａ和Ｂ是总的Ｈ和b
bool LinearAlignment(map<double, ImageFrame> &all_image_frame, Vector3d &g, VectorXd &x)//求尺度s，重力加速度g，速度v
{
    int all_frame_count = all_image_frame.size();
    int n_state = all_frame_count * 3 + 3 + 1;//这里的重力加速度是３维的

    MatrixXd A{n_state, n_state};
    A.setZero();
    VectorXd b{n_state};
    b.setZero();

    map<double, ImageFrame>::iterator frame_i;
    map<double, ImageFrame>::iterator frame_j;
    int i = 0;
    for (frame_i = all_image_frame.begin(); next(frame_i) != all_image_frame.end(); frame_i++, i++)
    {
        frame_j = next(frame_i);

        MatrixXd tmp_A(6, 10);//Ｈ矩阵的行列数
        tmp_A.setZero();//tmp_A和tmp_b是相邻间的速度的临时变量
        VectorXd tmp_b(6);
        tmp_b.setZero();

        double dt = frame_j->second.pre_integration->sum_dt;

        tmp_A.block<3, 3>(0, 0) = -dt * Matrix3d::Identity();
        tmp_A.block<3, 3>(0, 6) = frame_i->second.R.transpose() * dt * dt / 2 * Matrix3d::Identity();
        tmp_A.block<3, 1>(0, 9) = frame_i->second.R.transpose() * (frame_j->second.T - frame_i->second.T) / 100.0;     
        tmp_b.block<3, 1>(0, 0) = frame_j->second.pre_integration->delta_p + frame_i->second.R.transpose() * frame_j->second.R * TIC[0] - TIC[0];
        //cout << "delta_p   " << frame_j->second.pre_integration->delta_p.transpose() << endl;
        tmp_A.block<3, 3>(3, 0) = -Matrix3d::Identity();
        tmp_A.block<3, 3>(3, 3) = frame_i->second.R.transpose() * frame_j->second.R;
        tmp_A.block<3, 3>(3, 6) = frame_i->second.R.transpose() * dt * Matrix3d::Identity();
        tmp_b.block<3, 1>(3, 0) = frame_j->second.pre_integration->delta_v;
        //cout << "delta_v   " << frame_j->second.pre_integration->delta_v.transpose() << endl;

        Matrix<double, 6, 6> cov_inv = Matrix<double, 6, 6>::Zero();
        //cov.block<6, 6>(0, 0) = IMU_cov[i + 1];
        //MatrixXd cov_inv = cov.inverse();
        cov_inv.setIdentity();//直接将协方差设置成了单位矩阵

        MatrixXd r_A = tmp_A.transpose() * cov_inv * tmp_A;
        VectorXd r_b = tmp_A.transpose() * cov_inv * tmp_b;

        A.block<6, 6>(i * 3, i * 3) += r_A.topLeftCorner<6, 6>();//这几步不知道是什么意思
        b.segment<6>(i * 3) += r_b.head<6>();

        A.bottomRightCorner<4, 4>() += r_A.bottomRightCorner<4, 4>();
        b.tail<4>() += r_b.tail<4>();

        A.block<6, 4>(i * 3, n_state - 4) += r_A.topRightCorner<6, 4>();
        A.block<4, 6>(n_state - 4, i * 3) += r_A.bottomLeftCorner<4, 6>();
    }
    A = A * 1000.0;
    b = b * 1000.0;
    x = A.ldlt().solve(b);
    double s = x(n_state - 1) / 100.0;
    ROS_DEBUG("estimated scale: %f", s);
    g = x.segment<3>(n_state - 4);
    ROS_DEBUG_STREAM(" result g     " << g.norm() << " " << g.transpose());
    if(fabs(g.norm() - G.norm()) > 1.0 || s < 0)//重力加速度的模不能相差太大
    {
        return false;
    }

    RefineGravity(all_image_frame, g, x);//重新计算重力加速度方向，是２维的，得到最优解，
    s = (x.tail<1>())(0) / 100.0;//状态向量的倒数第一个是尺度因子
    (x.tail<1>())(0) = s;
    ROS_DEBUG_STREAM(" refine     " << g.norm() << " " << g.transpose());
    if(s < 0.0 )
        return false;   
    else
        return true;
}

bool VisualIMUAlignment(map<double, ImageFrame> &all_image_frame, Vector3d* Bgs, Vector3d &g, VectorXd &x)
{
    solveGyroscopeBias(all_image_frame, Bgs);

    if(LinearAlignment(all_image_frame, g, x))
        return true;
    else 
        return false;
}
