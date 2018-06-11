#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "json.hpp"
#include "spline.h"

using namespace std;

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_first_of("}");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

double distance(double x1, double y1, double x2, double y2)
{
	return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}
int ClosestWaypoint(double x, double y, const vector<double> &maps_x, const vector<double> &maps_y)
{

	double closestLen = 100000; //large number
	int closestWaypoint = 0;

	for(int i = 0; i < maps_x.size(); i++)
	{
		double map_x = maps_x[i];
		double map_y = maps_y[i];
		double dist = distance(x,y,map_x,map_y);
		if(dist < closestLen)
		{
			closestLen = dist;
			closestWaypoint = i;
		}

	}

	return closestWaypoint;

}

int NextWaypoint(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{

	int closestWaypoint = ClosestWaypoint(x,y,maps_x,maps_y);

	double map_x = maps_x[closestWaypoint];
	double map_y = maps_y[closestWaypoint];

	double heading = atan2((map_y-y),(map_x-x));

	double angle = fabs(theta-heading);
  angle = min(2*pi() - angle, angle);

  if(angle > pi()/4)
  {
    closestWaypoint++;
  if (closestWaypoint == maps_x.size())
  {
    closestWaypoint = 0;
  }
  }

  return closestWaypoint;
}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int next_wp = NextWaypoint(x,y, theta, maps_x,maps_y);

	int prev_wp;
	prev_wp = next_wp-1;
	if(next_wp == 0)
	{
		prev_wp  = maps_x.size()-1;
	}

	double n_x = maps_x[next_wp]-maps_x[prev_wp];
	double n_y = maps_y[next_wp]-maps_y[prev_wp];
	double x_x = x - maps_x[prev_wp];
	double x_y = y - maps_y[prev_wp];

	// find the projection of x onto n
	double proj_norm = (x_x*n_x+x_y*n_y)/(n_x*n_x+n_y*n_y);
	double proj_x = proj_norm*n_x;
	double proj_y = proj_norm*n_y;

	double frenet_d = distance(x_x,x_y,proj_x,proj_y);

	//see if d value is positive or negative by comparing it to a center point

	double center_x = 1000-maps_x[prev_wp];
	double center_y = 2000-maps_y[prev_wp];
	double centerToPos = distance(center_x,center_y,x_x,x_y);
	double centerToRef = distance(center_x,center_y,proj_x,proj_y);

	if(centerToPos <= centerToRef)
	{
		frenet_d *= -1;
	}

	// calculate s value
	double frenet_s = 0;
	for(int i = 0; i < prev_wp; i++)
	{
		frenet_s += distance(maps_x[i],maps_y[i],maps_x[i+1],maps_y[i+1]);
	}

	frenet_s += distance(0,0,proj_x,proj_y);

	return {frenet_s,frenet_d};

}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(double s, double d, const vector<double> &maps_s, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int prev_wp = -1;

	while(s > maps_s[prev_wp+1] && (prev_wp < (int)(maps_s.size()-1) ))
	{
		prev_wp++;
	}

	int wp2 = (prev_wp+1)%maps_x.size();

	double heading = atan2((maps_y[wp2]-maps_y[prev_wp]),(maps_x[wp2]-maps_x[prev_wp]));
	// the x,y,s along the segment
	double seg_s = (s-maps_s[prev_wp]);

	double seg_x = maps_x[prev_wp]+seg_s*cos(heading);
	double seg_y = maps_y[prev_wp]+seg_s*sin(heading);

	double perp_heading = heading-pi()/2;

	double x = seg_x + d*cos(perp_heading);
	double y = seg_y + d*sin(perp_heading);

	return {x,y};

}

vector<string> successor_states(string state, int lane) 
{
    /*
    Provides the possible next states given the current state for the FSM 
    discussed in the course, with the exception that lane changes happen 
    instantaneously, so LCL and LCR can only transition back to KL.
    */
	int lanes_available = 3;
  vector<string> states;
  states.push_back("KL");

  if(state.compare("KL") == 0) 
  {
    states.push_back("PLCL");
    states.push_back("PLCR");
  } 
  else if (state.compare("PLCL") == 0) 
  {
    if (lane > 0) 
    {
      states.push_back("PLCL");
      states.push_back("LCL");
    }
  } 
  else if (state.compare("PLCR") == 0) 
  {
    if (lane < lanes_available - 1) 
    {
      states.push_back("PLCR");
      states.push_back("LCR");
    }
  }
  //If state is "LCL" or "LCR", then just return "KL"
  return states;
}

vector<double> lane_speeds(vector<double> left_lane, vector<double> middle_lane, vector<double> right_lane)
{
	double left_lane_speed = 0.0;
	double middle_lane_speed = 0.0;
	double right_lane_speed = 0.0;
	vector<double> avg_lane_speeds = {0.0, 0.0, 0.0};
	
	for (int i = 0; i < left_lane.size(); i++)
	{
		left_lane_speed += left_lane[i];
	}
	avg_lane_speeds[0] = left_lane_speed / left_lane.size();
	
	for (int i = 0; i < middle_lane.size(); i++)
	{
		middle_lane_speed += middle_lane[i];
	}
	avg_lane_speeds[1] = middle_lane_speed / middle_lane.size();
	
	for (int i = 0; i < right_lane.size(); i++)
	{
		right_lane_speed += right_lane[i];
	}
	avg_lane_speeds[2] = right_lane_speed / right_lane.size();
	
	return avg_lane_speeds;
}

vector<double> lane_costs(vector<double> avg_lane_speed, vector<int> num_cars_in_lane)
{

	vector<double> lane_cost = {0.0, 0.0, 0.0};
	// higer cost for slower lanes
	for (int i = 0; i < avg_lane_speed.size(); i++)
	{
		lane_cost[i] += (1 / (avg_lane_speed[i]));
	}
	
	//Higer cost for crowded lanes
	for (int i = 0; i < num_cars_in_lane.size(); i++)
	{
		lane_cost[i] += (10 * num_cars_in_lane[i]);
	}
	
	return lane_cost;
}

string best_state(int lane, vector<double> lane_costs, vector<string> safe_states)
{
	//string best_next_state; //= None
    
	double min_cost = 9999999;
	double temp_min_cost = 0;
	int best_cost_lane;

	for(int i = 0; i < lane_costs.size(); i++)
	{
		//state = possible_successor_states[i]
    temp_min_cost  = lane_costs[i];
    if (temp_min_cost < min_cost)
		{
      min_cost = temp_min_cost;
			best_cost_lane = i;
		}
	}

	if (best_cost_lane == lane)
	{
		return "KL";
	}

	for(int i = 0; i < safe_states.size(); i++)
	{
		if(lane == 0) // left lane
		{
			if(best_cost_lane == 1 && safe_states[i].compare("LCR") == 0)
			{
				return "LCR";
			}
			if (best_cost_lane == 2)
			{
				if(lane_costs[1] < lane_costs[lane])
				{
					best_cost_lane = 1;
				}
			}
		}

		if(lane == 1) // middle lane
		{
			if(best_cost_lane == 0 && safe_states[i].compare("LCL") == 0)
			{
				return "LCL";
			}
			else if(best_cost_lane == 2 && safe_states[i].compare("LCR") == 0)
			{
				return "LCR";
			}
		}

		if(lane == 2) // right lane
		{
			if (best_cost_lane == 0)
			{
				if(lane_costs[1] < lane_costs[lane])
				{
					best_cost_lane = 1;
				}
			}
			if(best_cost_lane == 1 && safe_states[i].compare("LCL") == 0)
			{
				return "LCL";
			}	
		}
	}

  for(int i = 0; i < safe_states.size(); i++)
	{
		if(lane == 0) // left lane
		{
			if(best_cost_lane == 1 && safe_states[i].compare("PLCR") == 0)
			{
				return "PLCR";
			}
			if (best_cost_lane == 2)
			{
				if(lane_costs[1] < lane_costs[lane])
				{
					best_cost_lane = 1;
				}
			}
		}

		if(lane == 1) // middle lane
		{
			if(best_cost_lane == 0 && safe_states[i].compare("PLCL") == 0)
			{
				return "PLCL";
			}
			else if(best_cost_lane == 2 && safe_states[i].compare("PLCR") == 0)
			{
				return "PLCR";
			}
		}

		if(lane == 2) // right lane
		{
			if (best_cost_lane == 0)
			{
				if(lane_costs[1] < lane_costs[lane])
				{
					best_cost_lane = 1;
				}
			}
			if(best_cost_lane == 1 && safe_states[i].compare("PLCL") == 0)
			{
				return "PLCL";
			}
		}
	}
	return "KL";
}


int main() {
  uWS::Hub h;

  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  vector<double> map_waypoints_x;
  vector<double> map_waypoints_y;
  vector<double> map_waypoints_s;
  vector<double> map_waypoints_dx;
  vector<double> map_waypoints_dy;

  // Waypoint map to read from
  string map_file_ = "../data/highway_map.csv";
  // The max s value before wrapping around the track back to 0
  double max_s = 6945.554;

  ifstream in_map_(map_file_.c_str(), ifstream::in);

  string line;
  while (getline(in_map_, line)) {
  	istringstream iss(line);
  	double x;
  	double y;
  	float s;
  	float d_x;
  	float d_y;
  	iss >> x;
  	iss >> y;
  	iss >> s;
  	iss >> d_x;
  	iss >> d_y;
  	map_waypoints_x.push_back(x);
  	map_waypoints_y.push_back(y);
  	map_waypoints_s.push_back(s);
  	map_waypoints_dx.push_back(d_x);
  	map_waypoints_dy.push_back(d_y);
  }

	// Start in lane 1
	int lane =1;

	// Have a reference velopcity to target
	double ref_vel = 0.0; //mph

	vector<double> left_lane_close_objs, middle_lane_close_objs, right_lane_close_objs;
	vector<double> speed_left_lane, speed_middle_lane, speed_right_lane;
	string ego_state = "kL";

  h.onMessage([&map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy, &lane, &ref_vel, &left_lane_close_objs, &middle_lane_close_objs, &right_lane_close_objs, &speed_left_lane, &speed_middle_lane, &speed_right_lane, &ego_state](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    //auto sdata = string(data).substr(0, length);
    //cout << sdata << endl;
    if (length && length > 2 && data[0] == '4' && data[1] == '2') {

      auto s = hasData(data);

      if (s != "") {
        auto j = json::parse(s);
        
        string event = j[0].get<string>();
        
        if (event == "telemetry") {
          // j[1] is the data JSON object
          
        	// Main car's localization Data
          double car_x = j[1]["x"];
          double car_y = j[1]["y"];
          double car_s = j[1]["s"];
          double car_d = j[1]["d"];
          double car_yaw = j[1]["yaw"];
          double car_speed = j[1]["speed"];

          // Previous path data given to the Planner
          auto previous_path_x = j[1]["previous_path_x"];
          auto previous_path_y = j[1]["previous_path_y"];
          // Previous path's end s and d values 
          double end_path_s = j[1]["end_path_s"];
          double end_path_d = j[1]["end_path_d"];

          // Sensor Fusion Data, a list of all other cars on the same side of the road.
          auto sensor_fusion = j[1]["sensor_fusion"];

          json msgJson;

          vector<double> next_x_vals;
          vector<double> next_y_vals;
			    vector<double> avg_lane_speeds;


          // TODO: define a path made up of (x,y) points that the car will visit sequentially every .02 seconds
			    int prev_size = previous_path_x.size();

			if (prev_size > 0)
			{
				car_s = end_path_s;
			}

			bool too_close = false;
			//bool FCW = false;
			bool left_lane_free = false;
			bool right_lane_free = false;
			double CIPV_vel;

			// Find ref_v to use
			for (int i = 0; i < sensor_fusion.size(); i++)
			{
				// Car is in my lane
				int obj_id = sensor_fusion[i][0];
				double vx = sensor_fusion[i][3];
				double vy = sensor_fusion[i][4];
				double check_speed = sqrt(vx*vx + vy*vy);
				double check_car_s = sensor_fusion[i][5];
				float d = sensor_fusion[i][6];
				int obj_fused_lane;

				check_car_s += ((double)prev_size * 0.02 * check_speed); // if using previous points can project s value out
				
				//if(d < (2+4*lane+2) && d > (2+4*lane-2))
				if(d < 4 && d > 0)
				{
					obj_fused_lane = 0;
				}
				else if(d < 8 && d > 4)
				{
					obj_fused_lane = 1;
				}
				else if(d < 12 && d > 8)
				{
					obj_fused_lane = 2;
				}
					
				//check s values greater than mine and s gap
				//if((check_car_s > car_s) && ((check_car_s-car_s) < 30))
        if((check_car_s - car_s) < 25.0 && ((check_car_s-car_s) > -10.0))
				{
					if(obj_fused_lane == 0 && find(left_lane_close_objs.begin(), left_lane_close_objs.end(), obj_id) == left_lane_close_objs.end())
					{
						left_lane_close_objs.push_back(obj_id);
					}
					else if(obj_fused_lane == 1 && find(middle_lane_close_objs.begin(), middle_lane_close_objs.end(), obj_id) == middle_lane_close_objs.end())
					{
						middle_lane_close_objs.push_back(obj_id);
					}
					else if(obj_fused_lane == 2 && find(right_lane_close_objs.begin(), right_lane_close_objs.end(), obj_id) == right_lane_close_objs.end())
					{
						right_lane_close_objs.push_back(obj_id);
					}
					//ref_vel = 29.5; //mph
					//too_close = true;
					//if(lane > 0)
					//{
					//	lane = 0;
					//}
				}
				else
				{
					if(obj_fused_lane == 0)
					{
						left_lane_close_objs.erase(remove(left_lane_close_objs.begin(), left_lane_close_objs.end(), obj_id), left_lane_close_objs.end());
					}
					else if(obj_fused_lane == 1)
					{
						middle_lane_close_objs.erase(remove(middle_lane_close_objs.begin(), middle_lane_close_objs.end(), obj_id), middle_lane_close_objs.end());
					}
					else if(obj_fused_lane == 2)
					{
						right_lane_close_objs.erase(remove(right_lane_close_objs.begin(), right_lane_close_objs.end(), obj_id), right_lane_close_objs.end());
					}
				}

				//if((check_car_s > car_s) && ((check_car_s-car_s) < 100))
        if((check_car_s - car_s) < 100.0 && ((check_car_s-car_s) > -15.0))
				{
					if(obj_fused_lane == 0)
					{
						speed_left_lane.push_back(check_speed);
					}
					else if(obj_fused_lane == 1)
					{
						speed_middle_lane.push_back(check_speed);
					}
					else if(obj_fused_lane == 2)
					{
						speed_right_lane.push_back(check_speed);
					}
					
					avg_lane_speeds = lane_speeds(speed_left_lane, speed_middle_lane, speed_right_lane);
				}
				else
				{
					if(obj_fused_lane == 0)
					{
						speed_left_lane.erase(remove(speed_left_lane.begin(), speed_left_lane.end(), check_speed), speed_left_lane.end());
					}
					else if(obj_fused_lane == 1)
					{
						speed_middle_lane.erase(remove(speed_middle_lane.begin(), speed_middle_lane.end(), check_speed), speed_middle_lane.end());
					}
					else if(obj_fused_lane == 2)
					{
						speed_right_lane.erase(remove(speed_right_lane.begin(), speed_right_lane.end(), check_speed), speed_right_lane.end());
					}

					avg_lane_speeds = lane_speeds(speed_left_lane, speed_middle_lane, speed_right_lane);
				}

				if(left_lane_close_objs.size() == 0)
				{
					speed_left_lane = {50.0};
				}

				if(middle_lane_close_objs.size() == 0)
				{
					speed_middle_lane = {50.0};
				}

				if(right_lane_close_objs.size() == 0)
				{
					speed_right_lane = {50.0};
				}

				if(d < (2+4*lane+2) && d > (2+4*lane-2))
				{
					//check s values greater than mine and s gap
					if((check_car_s > car_s) && ((check_car_s-car_s) < 30))
					{
						too_close = true;
						CIPV_vel = check_speed;

						vector<string> possible_states = successor_states(ego_state, lane);
						vector<string> safe_states;
						string final_best_state;

						// Determine is adjacent lanes are occupied
						if(lane == 0 && middle_lane_close_objs.size() == 0)
						{
							right_lane_free = true;
						}
						else if(lane == 1)
						{
							if(left_lane_close_objs.size() == 0)
							{
								left_lane_free = true;
							}

							if(right_lane_close_objs.size() == 0)
							{
								right_lane_free = true;
							}
						}
						else if(lane == 2 && middle_lane_close_objs.size() == 0)
						{
							left_lane_free = true;
						}

						for(int i = 0; i < possible_states.size(); i++)
						{
							if(possible_states[i].compare("KL") == 0) 
							{
								safe_states.push_back(possible_states[i]);
							} 
							
							if (possible_states[i].compare("PLCL") == 0) 
							{
								if (left_lane_free) 
								{
									safe_states.push_back(possible_states[i]);
								}
							} 
							
							if (possible_states[i].compare("PLCR") == 0) {
								if (right_lane_free) 
								{
									safe_states.push_back(possible_states[i]);
								}
							}

							if (possible_states[i].compare("LCL") == 0) 
							{
								if (left_lane_free) 
								{
									safe_states.push_back(possible_states[i]);
								}
							} 
							
							if (possible_states[i].compare("LCR") == 0) {
								if (right_lane_free) 
								{
									safe_states.push_back(possible_states[i]);
								}
							}
						}

						vector<int> num_cars_in_lanes;
						num_cars_in_lanes.push_back(left_lane_close_objs.size());
						num_cars_in_lanes.push_back(middle_lane_close_objs.size());
						num_cars_in_lanes.push_back(right_lane_close_objs.size());

						vector<double> lanes_cost = lane_costs(avg_lane_speeds, num_cars_in_lanes);

            //vector<double> lanes_cost = get_lane_costs(lane, avg_lane_speeds, safe_states, num_cars_in_lanes); //Calculate Lane costs

						final_best_state = best_state(lane, lanes_cost, safe_states);

						if(final_best_state == "PLCL")
						{

						}
						else if(final_best_state == "LCL")
						{
							lane -= 1;
						}
						else if(final_best_state == "PLCR")
						{

						}
						else if(final_best_state == "LCR")
						{
							lane +=1;
						}
						else if(final_best_state == "KL")
						{
							
						}
						
						ego_state = final_best_state;
					}
				}
			}

			/*if(too_close && ref_vel/2.24 > CIPV_vel + 0.5)
			{
				ref_vel -= 0.224;
			}
			else if(too_close && ref_vel/2.24 >= CIPV_vel)
			{
				ref_vel -= 0.5;
			}*/
      if (too_close)
      {
        ref_vel -= 0.4;
      }
			else if(ref_vel < 49.5)
			{
				ref_vel += 0.224;
			}

			// Create a list of widely spaced (x,y) waypoints, evenly spaced at 30m
			// Later we will interoplate these waypoints with a spline and fill it in with more points that we can control speed
			vector<double> ptsx;
			vector<double> ptsy;

			// reference x,y,yaw states
			// either we will reference the starting points as where the car is or at the previous path and points
			double ref_x = car_x;
			double ref_y = car_y;
			double ref_yaw = deg2rad(car_yaw);

			// if previous size is almost empty, use the car as starting reference
			if(prev_size < 2)
			{
				// Use two points that make the path tangent to the car
				double prev_car_x = car_x - cos(car_yaw);
				double prev_car_y = car_y - sin(car_yaw);

				ptsx.push_back(prev_car_x);
				ptsx.push_back(car_x);

				ptsy.push_back(prev_car_y);
				ptsy.push_back(car_y);
			}
			// use the previouss path' end points as starting reference
			else
			{
				ref_x = previous_path_x[prev_size-1];
				ref_y = previous_path_y[prev_size-1];

				double ref_x_prev = previous_path_x[prev_size-2];
				double ref_y_prev = previous_path_y[prev_size-2];
				ref_yaw = atan2(ref_y - ref_y_prev, ref_x - ref_x_prev);

				// Use two points that make the path tangent to the previouse path's end poitns
				ptsx.push_back(ref_x_prev);
				ptsx.push_back(ref_x);

				ptsy.push_back(ref_y_prev);
				ptsy.push_back(ref_y);
			}

			// In Frenet add evenly 30m spaced points ahead of the starting reference
			vector<double> next_wp0 = getXY(car_s+30, (2+4*lane), map_waypoints_s, map_waypoints_x, map_waypoints_y);
			vector<double> next_wp1 = getXY(car_s+60, (2+4*lane), map_waypoints_s, map_waypoints_x, map_waypoints_y);
			vector<double> next_wp2 = getXY(car_s+90, (2+4*lane), map_waypoints_s, map_waypoints_x, map_waypoints_y);

			ptsx.push_back(next_wp0[0]);
			ptsx.push_back(next_wp1[0]);
			ptsx.push_back(next_wp2[0]);

			ptsy.push_back(next_wp0[1]);
			ptsy.push_back(next_wp1[1]);
			ptsy.push_back(next_wp2[1]);

			for (int i = 0; i < ptsx.size(); i++)
			{
				// Shift car reference angle to 0 degrees
				double shift_x = ptsx[i]-ref_x;
				double shift_y = ptsy[i]-ref_y;

				ptsx[i] = (shift_x * cos(0 - ref_yaw) - shift_y * sin(0 - ref_yaw));
				ptsy[i] = (shift_x * sin(0 - ref_yaw) + shift_y * cos(0 - ref_yaw));
			}

			// Create a Spline
			tk::spline s;

			// Set (x,y) points to spline
			s.set_points(ptsx,ptsy);

			// Define the actual (x,y) points we will use or the planner
			//vector<double> next_x_vals;
			//vector<double> next_y_vals;

			// Start with all of the previous path points from last time
			for (int i = 0; i < previous_path_x.size(); i++)
			{
				next_x_vals.push_back(previous_path_x[i]);
				next_y_vals.push_back(previous_path_y[i]);
			}

			// Calculate how to break up Spline points so that we travel at our desired reference velocity
			double target_x = 30.0;
			double target_y = s(target_x);
			double target_dist = sqrt((target_x) * (target_x) + (target_y) * (target_y));
						
			double x_add_on = 0;

			// Fill up the rest of our path planner after filling it with previous points, here we will always putput 50 points
			for (int i = 1; i <= 50-previous_path_x.size(); i++)
			{
				double N = (target_dist/(0.02*ref_vel/2.24));
				double x_point = x_add_on + (target_x) / N;
				double y_point = s(x_point);

				x_add_on = x_point;

				double x_ref = x_point;
				double y_ref = y_point;

				// rotate back to normal after rotating it eariler
				x_point = (x_ref * cos(ref_yaw) - y_ref * sin(ref_yaw));
				y_point = (x_ref * sin(ref_yaw) + y_ref * cos(ref_yaw));

				x_point += ref_x;
				y_point += ref_y;

				next_x_vals.push_back(x_point);
				next_y_vals.push_back(y_point);
			}

	/*
							double dist_inc = 0.3;
							for(int i = 0; i < 50; i++)
							{
									double next_s = car_s+(i+1)*dist_inc;
									double next_d = 6;
									vector<double> xy = getXY(next_s, next_d, map_waypoints_s, map_waypoints_x, map_waypoints_y);
									
									next_x_vals.push_back(xy[0]);
									next_y_vals.push_back(xy[1]);
							}*/

			msgJson["next_x"] = next_x_vals;
			msgJson["next_y"] = next_y_vals;

			auto msg = "42[\"control\","+ msgJson.dump()+"]";

			//this_thread::sleep_for(chrono::milliseconds(1000));
			ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
			
        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}
