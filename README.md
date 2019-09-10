# Image-Alignment-using-Inertial-Sensor-UKF

Merging short-exposure frames can provide an image with reduced noise in low light conditions. However, how best to align images before merging is an open problem. To improve the performance of alignment, we propose an inertia-sensor aided algorithm for smartphone burst photography, which takes rotation and out-plane relative movement into account. To calculate homography between frames, a three by three rotation matrix is calculated from gyro data recorded by smartphone inertia sensor and three dimensional translation vector are estimated by matched feature points detected from two frames. The rotation matrix and translations are com- bined to form the initial guess of homography. An unscented Kalman filter is utilized to provide a more accurate homogra- phy estimation. We have tested the algorithm on a variety of different scenes with different camera relative motions. We compare the proposed method to benchmark single-image and multi-image denoising methods with favourable results.

## How to use
1. Put a sequence of noisy **RAW** frames and its corresponding inertial sensor data in the current folder
2. Run **Main_lowLight.m** to get merged (denoised) frame in RAW format
3. Run **ProcessMergedImage.m** to get PNG result.

## Citation
@article{InertiaUKF:2018,

  title={Inertia Sensor Aided Alignment for Burst Pipeline in Low Light Conditions},
  
  author={Shuang Zhang and Robert L. Stevenson},
  
  journal={2018 25th IEEE International Conference on Image Processing (ICIP)},
  
  year={2018},
  
  pages={3953-3957}
}
