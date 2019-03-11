/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h>
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

#define EPS 0.00001

std::default_random_engine gen;

using namespace std;
using std::string;
using std::vector;
using std::normal_distribution;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
    std::default_random_engine gen;

    num_particles = 100;  // TODO: Set the number of particles
    double std_x, std_y, std_theta;
    std_x=std[0];
    std_y=std[1];
    std_theta=std[2];
    normal_distribution<double> noise_x(x,std_x);
    normal_distribution<double> noise_y(y,std_y);
    normal_distribution<double> noise_theta(theta,std_theta);
    

    for (int i=0;i<num_particles;++i){
        
        // create the initial x, y, theta with normalized gaussian distribution
        double x_ini, y_ini, theta_ini;
        x_ini=noise_x(gen);
        y_ini=noise_y(gen);
        theta_ini=noise_theta(gen);
        Particle particle;
        particle.id=i;
        particle.x=x_ini;
        particle.y=y_ini;
        particle.theta=theta_ini;
        particle.weight=1.0;
        
        particles.push_back(particle);
        //weights.push_back(1.0);
    }
    is_initialized=true;
 


}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */
    std::default_random_engine gen;
    double x_pre, y_pre, theta_pre;
    //define the nosie of each state
    normal_distribution<double> dist_x(0,std_pos[0]);
    normal_distribution<double> dist_y(0,std_pos[1]);
    normal_distribution<double> d(0,std_pos[2]);
    for(int i=0;i<num_particles;++i){
        //define gaussian noise
        x_pre=particles[i].x;
        y_pre=particles[i].y;
        theta_pre=particles[i].theta;
        
        if(abs(yaw_rate)>0.001){
            particles[i].x=x_pre+velocity/yaw_rate*(sin(theta_pre+yaw_rate*delta_t)-sin(theta_pre));
            particles[i].y=y_pre+velocity/yaw_rate*(cos(theta_pre)-cos(theta_pre+yaw_rate*delta_t));
            particles[i].theta=theta_pre+yaw_rate*delta_t;
        }
        else{
            particles[i].x=x_pre+velocity*cos(theta_pre)*delta_t;
            particles[i].y=y_pre+velocity*sin(theta_pre)*delta_t;
            //particles[i].theta=theta_pre;
        }
        //Adding random noise to the state
        
        particles[i].x+=dist_x(gen);
        particles[i].y+=dist_y(gen);
        particles[i].theta+=d(gen);
    }
    

}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */
    // leave the transformation to the next function
    for(unsigned int i=0;i<observations.size();++i){
        double curr_dist;
        double min_dist=1000;
        int closet_id;
        double x_mark, y_mark;
        for(unsigned int j=0;j<predicted.size();++j){
            curr_dist=dist(predicted[j].x, predicted[j].y,observations[i].x, observations[i].y);
            if(curr_dist<min_dist){
                closet_id=predicted[j].id;
                min_dist=curr_dist;
                x_mark=predicted[j].x;
                y_mark=predicted[j].y;
            }
        
        }
        observations[i].id=closet_id;
        //observations[i].x=x_mark;
        //observations[i].y=y_mark;
    }

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                    vector<LandmarkObs> observations, 
                                    Map map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */
   
    //Predict measurements to all map landmarks within sensor range
    for (int i=0;i<num_particles;++i){
        vector<LandmarkObs> obs_landmark; //created a vector of observed meas with all map_landmarks
        //vector<LandmarkObs> obs_par;
        particles[i].weight=1.0;
        for(unsigned int m=0;m<map_landmarks.landmark_list.size();++m){
            double curr_dist;
            int within_range;
          double par_x=particles[i].x;
          double par_y=particles[i].y;
          double Lm_x=map_landmarks.landmark_list[m].x_f;
          double Lm_y=map_landmarks.landmark_list[m].y_f;
            curr_dist=dist(Lm_x, Lm_y,par_x, par_y);
            if(curr_dist<=sensor_range){
              within_range=map_landmarks.landmark_list[m].id_i;
              LandmarkObs obs={within_range,Lm_x, Lm_y};
              obs_landmark.push_back(obs);
                
            }
        }
        vector<LandmarkObs> trans_obs;
        //vector<LandmarkObs> obs_global; //create another vector to store transformed observations
        // transform obs to global coordinates
        double map_x, map_y, theta_m;
        theta_m=particles[i].theta;
        for(unsigned int j=0;j<observations.size();++j){
            map_x=cos(theta_m)*observations[j].x-sin(theta_m)*observations[j].y+particles[i].x;
    		map_y=sin(theta_m)*observations[j].x+cos(theta_m)*observations[j].y+particles[i].y;
    		LandmarkObs trans_pos;
    		trans_pos={observations[j].id, map_x, map_y};
            //trans_obs[j]=transform(particles[i],observations[j]);
            //obs_global.push_back(trans_pos);
            trans_obs.push_back(trans_pos);
        }
        //data association to store respective landmark observations
        dataAssociation(obs_landmark, trans_obs);
        double weight_par=1.0;
        //calculate weights
        int land_id;
        particles[i].weight=1.0;
        for(unsigned int o=0;o<trans_obs.size();++o){
            land_id=trans_obs[o].id;
            //Map::single_landmark_s landmark = map_landmarks.landmark_list.at(land_id-1);
            //double mu_x=landmark.x_f;
            //double mu_y=landmark.y_f;
            double mu_x, mu_y;
            for (unsigned int k = 0; k < obs_landmark.size(); k++) {
        		if (obs_landmark[k].id == land_id) {
          			mu_x = obs_landmark[k].x;
          			mu_y = obs_landmark[k].y;
        		}
      		}
            //weight_par=multi_gussan(std_landmark[0], std_landmark[1], trans_obs[n].x, trans_obs[n].y,
            //                         mu_x, mu_y);
            double x_term = pow(trans_obs[o].x - mu_x, 2) / (2 * pow(std_landmark[0], 2));
            double y_term = pow(trans_obs[o].y - mu_y, 2) / (2 * pow(std_landmark[1], 2));
            double prob= exp(-(x_term + y_term)) / (2 * M_PI * std_landmark[0] * std_landmark[1]);
            weight_par*=prob;
        }
        particles[i].weight=weight_par;
        //weights.push_back(weight_par);
        weights.push_back(particles[i].weight);
    }
    
    

}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
  std::random_device rd;
  std::mt19937 gen(rd());
  std::discrete_distribution<int> dist(weights.begin(), weights.end());
    vector<Particle> resample_par(particles.size());
    resample_par.resize(num_particles);
    for(int j=0;j<num_particles;++j){
        int index=dist(gen);
        resample_par[j]=particles[index];
    }
    particles=resample_par;
  weights.clear();
}

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
{
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

	//Clear the previous associations
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

	particle.associations= associations;
 	particle.sense_x = sense_x;
 	particle.sense_y = sense_y;

 	return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}