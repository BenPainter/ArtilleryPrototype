// ArtilleryPrototype.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
#define _USE_MATH_DEFINES
#include <iostream>
#include <math.h>   // sin,cos,sqrt,atan
#include <cmath>    // for Pi
#include <vector>
#include <cassert>
using namespace std;


// c - drag coeficient
// p - density of gas
// v - velocity
// a - surface area
double dragForce(const double c, const double p, const double v, const double a)
{
   // returns d (force in newtons)
   return (.5 * c * p * (v * v)  * a  );
}

double verticalComp(const double s, const double a)
{
   // returns dy
   return s * cos(a);
}

double horizontalComp(const double s, const double a)
{
   // returns dx
   return s * sin(a);
}


// v = r  c = d
double linearInter(const double d, const double d0, const double  r0, const double d1, const double r1)
{
   // returns c
   return r0 + ((d - d0) * (r1 - r0)) / (d1 - d0);

}

double speed(const double dx, const double dy)
{
   return sqrt((dx * dx) + (dy * dy));
}

// v - old velocity
// a - acceleration
// t - time
double computeVelocity(const double v, const double a, const double t)
{
   return v + (a * t);
}

// s - orignial position
// v - velocity
// a - acceleration
// t - time
double computeDistance(double s, double v, double a, double t)
{
   // returns the new position
   return s + (v * t) + (.5 * a * (t * t));
}

double computeAngle(double dx, double dy)
{
   double newAngle =  atan2(dx, dy);
   return newAngle; // -floor(newAngle / (2.00 * M_PI)) * (2.00 * M_PI);
}



double tableLookUp(const vector<double> chart1, const vector <double> chart2, const double d)
{
   assert(chart1.size() == chart2.size()); // the two charts should always be the same size

   int i = 0;
   assert(chart1[i] <= d);
   while (i < chart1.size())
   {
      if (chart1[i] == d)
      {
         return chart2[i];
      }
      else if (chart1[i] > d)
      {
         return linearInter(d, chart1[i - 1], chart2[i - 1], chart1[i], chart2[i]);
      }
      i++;
   }

   return 0.0;
}



int main()
{

   vector <double> altitudeArray = { 0.0, 1000.0, 2000.0, 3000.0, 4000.0, 5000.0, 6000.0, 7000.0, 8000.0, 9000.0, 10000.0, 15000.0, 20000.0, 25000.0 };
   vector <double> densityArray = {1.225, 1.112, 1.007, 0.9093, 0.8194, 0.7364, 0.6601, 0.59, 0.5258, 0.4671, .04135, 0.1948, 0.08891, 0.04008};
   vector <double> soundArray = {340.0, 336.0, 332.0, 328.0, 324.0, 320.0, 316.0, 312.0, 308.0, 303.0, 299.0, 295.0, 295.0, 295.0};
   vector <double> gravityArray = {9.807, 9.804, 9.801, 9.797, 9.794, 9.791, 9.788, 9.785, 9.782, 9.779, 9.776, 9.761, 9.745, 9.730};
   vector <double> machArray = {0.3, 0.5, 0.7, 0.89, 0.92, 0.96, 0.98, 1.0, 1.02, 1.06, 1.24, 1.53, 1.99, 2.87, 2.89, 5.0};
   vector <double> dragArray = {0.1629, 0.1659, 0.2031, 0.2597, 0.301, 0.3287, 0.4002, 0.4258, 0.4335, 0.4483, 0.4064, 0.3663, 0.2897, 0.2297, 0.2306, 0.2656};
   double x = 0.0;
   double y = 0.0;
   double t = 1;            // one second
   double angle = 0.5235;   // 10 degrees 0.174533
   double v = 827.0;        // m/s
   double mass = 46.7;      // kg
   double r = 0.077445;     // meters
   double area = 0.018842;  // meters2
   double p;
   double vSound;
   double vMach;
   double c;
   double force;
   double acc;
   double dx = horizontalComp(v, angle);
   double dy = verticalComp(v, angle);
   double ddy = 0.0;
   double ddx = 0.0;

   int i = 0;
   while (i <= 100)
   {

      p = 1.225;// = tableLookUp(altitudeArray, densityArray, y);      // air density
      vSound = 340;// tableLookUp(altitudeArray, soundArray, y);   // speed of sound
      vMach = v / vSound;                                   // Mach 
      c = tableLookUp(machArray, dragArray, vMach);         //linearInter(2.43235, 1.990, 0.2897, 2.870, 0.2297);
      force = dragForce(c, p, v, area);
      acc = (force / mass);
      ddx = horizontalComp(-acc, angle);
      ddy = -tableLookUp(altitudeArray, gravityArray, y) + verticalComp(-acc, angle); // gravity = -9.807
      cout << "x: " << x << endl;
     // cout << "y: " << y << endl;
      cout << "dx: " << dx << endl;
      //cout << "dy: " << dy << endl;
      cout << "ddx: " << ddx << endl;
      //cout << "ddy: " << ddy << endl;
      x = computeDistance(x, dx, ddx, t);
      y = computeDistance(y, dy, ddy, t);
      dx = computeVelocity(dx, ddx, t);
      dy = computeVelocity(dy, ddy, t);
      angle = computeAngle(dx, dy);
     // cout << force << endl;
      //cout << acc << endl;
     /* cout << "x: " << x << endl;
      cout << "y: " << y << endl;
      
      */
      cout << "angle: " << angle << endl;
      cout << "--------------------------\n";
      i++;
   }
   
}
