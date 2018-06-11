# CarND-Path-Planning-Project
Self-Driving Car Engineer Nanodegree Program
   
---

## Introduction
The goal of this project is to safely navigate around a virtual highway with other traffic that is driving +-10 MPH of the 50 MPH speed limit. The car's localization and sensor fusion data are provided, there is also a sparse map list of waypoints around the highway. The car should try to go as close as possible to the 50 MPH speed limit, which means passing slower traffic when possible, note that other cars will try to change lanes too. The car should avoid hitting other cars at all cost as well as driving inside of the marked road lanes at all times, unless going from one lane to another. The car should be able to make one complete loop around the 6946m highway. Since the car is trying to go 50 MPH, it should take a little over 5 minutes to complete 1 loop. Also the car should not experience total acceleration over 10 m/s^2 and jerk that is greater than 10 m/s^3.

## Implementation
To achieve the goal of this project, the implementation had to be broken down into different categories.

### 1. Prediction
In order for a vehicle to drive safely around a track, it needs to know where are the objects around it along with its location. I used and analyzed the localization and sensor fusion data for all vehicles to understand my surrounding including the ego vehicle. I then created a few flags `too_cloose`, `laft_lane_free`, `right_lane_free`, and `CIPV_vehicle`. Based on the sensor fusion and localization data, those flags helped me making trajectories decision since it is keep track if there is a vehicle too close to the ego vehicle, which vehicle is closest in my path, and if the left and right lane are free for a lane change. The code for this portion is included in Main.cpp from line 473 to line 660

### 2. Behavior
Now we know the ego vehicle's surroundings, it is time to plan the ego vehicle's path. I created a cost function that keeps track of the cost of all three lanes at all times. This cost function penalizes crowded lanes and lanes with slower traffic. based on the cost function, the planner uses a Finite State Machine to decide if a lane change is necessary and if so which lane, and also decide if the ego vehicle needs to accelerate or decelerate. The code for this portion is included in Main.cpp from line 662 - line 714

### 2. Trajectory
After a desired lane of the ego vehicle has been set, We create (x, y) coordinate of waypoints ahead of the ego vehicle using the corresponding spaced (s,d) Frenet coordinates. From these waypoints, we use a spline function to generate evenly spaced points for a smooth trajectory towards a desired horizon. The car is expected to visit each of these points every 0.02 secs, so the number of the points are calculated by `N = (target_dist/(0.02*ref_vel/2.24))` in order to travel at the desired reference velocity with a smooth trajectory. The code for this portion is included in Main.cpp from line 716 to line 826.

## Conclusion
Overall, I am very satisfied with the path planner performance. The ego vehicle was able to drive more than 43 miles around the track without any incidents. However, I think there is still room for improvements. One thing I would like to add in the future is a better predication algorithm, that is able to better predict target vehicles changing lanes so we can avoid another vehicle crashing into the ego vehicle.

[//]: # (Image References)

[image1]: ./screenshot.png "Track Screenshot.png"

![alt text][image1]

