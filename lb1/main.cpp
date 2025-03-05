#include <stdio.h>
#include <SDL2/SDL.h>
#include <iostream>
#include <cmath>
#include <algorithm> 
#include <vector>
#include <SDL2/SDL_ttf.h>
#include <string>
#include <iomanip>
#include <sstream>
#include <limits>
#include <fstream>
using namespace std;

class Vector{
public:
	double x;
	double y;
	double z;

	Vector(double x, double y, double z): x(x), y(y), z(z) {};

	Vector operator+(const Vector& vec) const{
		return Vector(x + vec.x, y + vec.y, z + vec.z);
	}

	Vector operator+(double value) const{
		return Vector(x + value, y + value, z + value);
	}

	Vector& operator+=(const Vector& vec){
		*this = *this + vec;
		return *this;
	}


	Vector operator*(double value) const{
		return Vector(x * value, y * value, z * value);
	}

	Vector operator/(double value) const{
		return Vector(x / value, y/ value, z / value);
	}

	Vector operator/(int value) const{
		return Vector(x / value, y/ value, z / value);
	}

	Vector& operator/=(int value){
		*this = *this / value;
		return *this;
	}

	Vector operator-(const Vector& vec) const{
		return Vector(x - vec.x, y - vec.y, z - vec.z);
	}

	Vector& operator-=(const Vector& vec){
		*this = *this - vec;
		return *this;
	}

	Vector operator-(double value) const{
		return Vector(x - value, y - value, z - value);
	}

	bool operator!=(const Vector& vec){
		if (x != vec.x || y != vec.y || z != vec.z) return true;
		return false;
	}
};

enum Method{Euler, MiddlePoint, RK_4, Trapezoid, PredictorCorrector, DOPRI_8};

double euclideanDistance(Vector& first, Vector& second){
	return sqrt( pow(first.x - second.x, 2) + pow(first.y - second.y, 2) + pow(first.z - second.z, 2) );
}


string methodToStr(Method method){
	switch (method){
		case Euler: return "Euler";
		case MiddlePoint: return "MiddlePoint";
		case RK_4: return "RK4";
		case Trapezoid: return "Trapezoid";
		case PredictorCorrector: return "PredictorCorrector";
		case DOPRI_8: return "DOPRI8";
	}
	return 0;
}

class Attractor{
public:
	int count_point;
	double b;
	vector <Vector> coordinates;
	double h;
	enum Orient {X, Y, Z};
	double lim = 0.001;

	struct Rotate{
		Orient orient;
		int angle;
	};

	Method current_method = PredictorCorrector;

	vector <Rotate> rotates;

	int height = 1000;
    int width = 1000;
    SDL_Window* window;
	SDL_Texture* textTexture;
	SDL_Surface* textSurface;
    SDL_Renderer* renderer;
	TTF_Font* font;
	
	double tol = 0.0001;

	Attractor (int count_point, double b, double h): count_point(count_point), b(b), h(h) {
		setCoordinates(&Attractor::predictorCorrectorMethod);
        
		SDL_Init(SDL_INIT_VIDEO);
        window = SDL_CreateWindow("Thomas Attractor",
            SDL_WINDOWPOS_UNDEFINED,
            SDL_WINDOWPOS_UNDEFINED,
            width, height,
            SDL_WINDOW_SHOWN);
        renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED | SDL_RENDERER_PRESENTVSYNC);
        SDL_SetRenderDrawBlendMode(renderer, SDL_BLENDMODE_BLEND);
        SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
        SDL_RenderClear(renderer);
        SDL_RenderPresent(renderer);

		if (TTF_Init() == -1) {
			std::cerr << "TTF_Init failed: " << TTF_GetError() << std::endl;
			return;
		}
		font = TTF_OpenFont("./Roboto-VariableFont_wdth,wght.ttf", 24);
		if (font == nullptr) {
			std::cerr << "Failed to load font: " << TTF_GetError() << std::endl;
			return;
		}
	};

	Vector d_dt(Vector point){
		double x = sin(point.y) - b * point.x;
		double y = sin(point.z) - b * point.y;
		double z = sin(point.x) - b * point.z;
		return Vector(x,y,z);
	}

	Vector findCenter(){
		Vector center(0,0,0);
		for (auto& point: coordinates){
			center += point;
		}
		center /= coordinates.size();
		return center;
	}

	void doSetCoordinates(){
		switch(current_method){
			case Euler:
				setCoordinates(&Attractor::EulerMethod);
				break;
			case MiddlePoint:
				setCoordinates(&Attractor::middlePointMethod);
				break;
			case RK_4:
				setCoordinates(&Attractor::RK4);
				break;
			case Trapezoid:
				setCoordinates(&Attractor::reverseTrapezoidMethod);
				break;			
			case PredictorCorrector:
				setCoordinates(&Attractor::predictorCorrectorMethod);
				break;
			case DOPRI_8:
				setCoordinates(&Attractor::DOPRI8);
				break;
			default:
				break;
		}
	}

	void setCoordinates(Vector (Attractor::*func)(Vector)){	
		coordinates.clear();
		coordinates.push_back(Vector(0.1, 0.2, 0.2));
		for (int i = 0; i < count_point; i++){
			if (i%10 == 0){
				coordinates.push_back((this->*func)(coordinates.back()));
			}
		}
	}

	void jacobian(double (&J)[3][3], Vector point){
		J[0][0] = 1 + (h*b/2);
		J[0][1] = -((h/2) * cos(point.y));
		J[0][2] = 0;
		J[1][0] = 0;
		J[1][1] = 1 + (h*b/2);
		J[1][2] = - ((h/2) * cos(point.z));
		J[2][0] = - ((h/2) * cos(point.x));
		J[2][1] = 0;
		J[2][2] = 1 + (h*b/2);
	}

	double determinate(double (&J)[3][3]) const{
		double det = J[0][0]*J[1][1]*J[2][2] + J[0][1]*J[1][2]*J[2][0] +
					 J[0][2]*J[1][0]*J[2][1] - J[0][2]*J[1][1]*J[2][0] -
					 J[0][1]*J[1][0]*J[2][2] - J[0][0]*J[1][2]*J[2][1];
		return det;
	}

	void inversionJacobian(double (&J)[3][3], double (&J_inv)[3][3]){
		double det = determinate(J);
		J_inv[0][0] = (J[1][1]*J[2][2] - J[1][2]*J[2][1]) / det;
		J_inv[0][1] = - ((J[1][0]*J[2][2] - J[2][0]*J[1][2]) / det);
		J_inv[0][2] = (J[1][0]*J[2][1] - J[2][0]*J[1][1]) / det;
		J_inv[1][0] = - ((J[0][1]*J[2][2] - J[2][1]*J[0][2]) / det);
		J_inv[1][1] = (J[0][0]*J[2][2] - J[0][2]*J[2][0]) / det;
		J_inv[1][2] = - ((J[0][0]*J[2][1] - J[0][1]*J[2][0]) / det);
		J_inv[2][0] = (J[0][1]*J[1][2] - J[0][2]*J[1][1]) / det;
		J_inv[2][1] = - ((J[0][0]*J[1][2] - J[1][0]*J[0][2]) / det);
		J_inv[2][2] = (J[0][0]*J[1][1] - J[0][1]*J[1][0]) / det;
	}

	//y(n+1) = y(n+1) - J^(-1)*F 
	void NewtonMethod(double (&J_inv)[3][3], Vector &new_point, Vector F){
		new_point.x = new_point.x - (J_inv[0][0]*F.x + J_inv[0][1]*F.y + J_inv[0][2]*F.z);
		new_point.y = new_point.y - (J_inv[1][0]*F.x + J_inv[1][1]*F.y + J_inv[1][2]*F.z);
		new_point.z = new_point.z - (J_inv[2][0]*F.x + J_inv[2][1]*F.y + J_inv[2][2]*F.z); 
	}

	Vector predictorCorrectorMethod(Vector point){
		Vector predict_point = EulerMethod(point);
		return doReverseTrapezoidMethod(point, predict_point);
	}

	Vector reverseTrapezoidMethod(Vector point){
		return doReverseTrapezoidMethod(point);
	}

	Vector doReverseTrapezoidMethod(Vector point, Vector predict_point = {0,0,0}){
		int max_iter = 100;

		//new_point == y(n+1)
		//point == y(n)
		Vector new_point = point;
		if (predict_point != point){
			new_point = predict_point;
		}

		Vector k1 = d_dt(point);
		for (int i = 0; i < max_iter; i++){
			double J[3][3], J_inv[3][3];

			Vector k2 = d_dt(new_point);
			Vector F = new_point - point + ((k1 + k2) * (h/2));

			jacobian(J, new_point);
			inversionJacobian(J, J_inv);

			NewtonMethod(J_inv, new_point, F);

			if (fabs(F.x) < lim && fabs(F.y) < lim && fabs(F.z) < lim) {
                break;
            }
		}
		return new_point;
	}

	Vector EulerMethod(Vector point){
		return point + d_dt(point) * h;
	}

	Vector middlePointMethod(Vector point){
		Vector k1 = d_dt(point);
		Vector k2 = d_dt(point + k1 * (h/2));
		return point + k2 * h;
	}

	Vector RK4(Vector point){
		Vector k1 = d_dt(point);
		Vector k2 = d_dt(point + k1 * (h/2));
		Vector k3 = d_dt(point + k2 * (h/2));
		Vector k4 = d_dt(point + k3 * h);
		return point + (k1/6 + k2/3 + k3/3 + k4/6) * h;
	}

	Vector DOPRI8(Vector point){
		Vector k1 = d_dt(point);
		Vector k2 = d_dt(point + k1*(h/5.0));
		Vector k3 = d_dt(point + k1*(h*3.0/40.0) + k2*(h*9.0/40.0));
		Vector k4 = d_dt(point + k1*(h*44.0/45.0) - k2*(h*56.0/15.0) + k3*(h*32.0/9.0));
		Vector k5 = d_dt(point + k1*(h*19372.0/6561.0) - k2*(h*25360.0/2187.0) + k3*(h*64448.0/6561.0) - k4*(h*212.0/729.0));
		Vector k6 = d_dt(point + k1*(h*9017.0/3168.0) - k2*(h*355.0/33.0) + k3*(h*46732.0/5247.0) + k4*(h*49.0/176.0) - k5*(h*5103.0/18656.0));
		Vector k7 = d_dt(point + k1*(h*35.0/384.0) + k3*(h*500.0/1113.0) + k4*(h*125.0/192.0) - k5*(h*2187.0/6784.0) + k6*(h*11.0/84.0));
		Vector add_point5 = point + k1*(h*35.0/384.0) + k3*(h*500.0/1113.0) + k4*(h*125.0/192.0) - k5*(h*2187.0/6784.0) + k6*(h*11.0/84.0);
		Vector add_point4 = point + k1*(h*5179.0/57600.0) + k3*(h*7571.0/16695.0) + k4*(h*393.0/640.0) - k5*(h*92097.0/339200.0) + k6*(h*187.0/2100.0) + k7*(h/40.0);
		
		Vector error_vector = add_point5 - add_point4;

		Vector null_vec = {0,0,0};
		double error = euclideanDistance(error_vector, null_vec);

		if (error < tol){
			h *= pow(tol/error, 1.0/9.0);
			return add_point5;
		} else{
			h *= pow (tol/error, 1.0/8.0);
			return DOPRI8(point);
		}

	}

	void drawPoint(size_t size, int x, int y){
        for(int i = 0; i != size; i++){
            for(int j = 0; j != size; j++){
                SDL_RenderDrawPoint(renderer, x+j, y+i);
            }
        }
    }

	void transformPoint(Vector& point, bool is_first, double k, double last_k = 1){
		if (is_first){
			point = point * k + (width / 2);
		}
		else{
			point = ((point - (width/2)) / last_k) * k + (width / 2);
		}
	}

	void transformAttractor(double k, bool is_first, double last_k = 1){
		for (auto& point: coordinates){
			transformPoint(point, is_first, k, last_k);
		}
	}

	void rotateAttractor(Rotate rotate){
		Vector center = findCenter();
		for (auto& point: coordinates){
			point -= center;
			switch(rotate.orient){
				case X:
					xRotate(point, rotate.angle);
					break;
				case Y:
					yRotate(point, rotate.angle);
					break;
				case Z:
					zRotate(point, rotate.angle);
					break;
				default:
					break;
			}
			point += center;
		}
	}

	void xRotate(Vector& point, double angle){
		double rad_angle = (angle * M_PI) / 180;
		double current_y;
		double current_z;

		current_y = point.y * cos(rad_angle) - point.z * sin(rad_angle);
		current_z = point.y * sin(rad_angle) + point.z * cos(rad_angle);

		point.z = current_z;
		point.y = current_y;
	}
	
	void yRotate(Vector& point, double angle){
		double rad_angle = (angle * M_PI) / 180;
		double current_z;
		double current_x;

		current_x = point.x * cos(rad_angle) + point.z * sin(rad_angle);
		current_z = point.x * (-sin(rad_angle)) + point.z * cos(rad_angle);

		point.x = current_x;
		point.z = current_z;
	}

	void zRotate(Vector& point, double angle){
		double rad_angle = (angle * M_PI) / 180;
		double current_y;
		double current_x;

		current_x = point.x * cos(rad_angle) - point.y * sin(rad_angle);
		current_y = point.x * sin(rad_angle) + point.y * cos(rad_angle);

		point.x = current_x;
		point.y = current_y;
	}

	void drawText(){
		SDL_Color textColor = {255, 255, 255, 255};

		stringstream ss;
		ss << "Method: " << methodToStr(current_method);
		textSurface = TTF_RenderText_Solid(font, ss.str().c_str(), textColor);
		textTexture = SDL_CreateTextureFromSurface(renderer, textSurface);

		SDL_FreeSurface(textSurface);
		SDL_Rect textRectM = {0, 0, textSurface->w, textSurface->h};
		SDL_RenderCopy(renderer, textTexture, nullptr, &textRectM);

		ss.str("");
		ss.clear();
		ss  << "b = " << std::fixed << std::setprecision(2) << b;
		textSurface = TTF_RenderText_Solid(font, ss.str().c_str(), textColor);
		textTexture = SDL_CreateTextureFromSurface(renderer, textSurface);

		SDL_FreeSurface(textSurface);
		SDL_Rect textRectb = {0, 30, textSurface->w, textSurface->h};
		SDL_RenderCopy(renderer, textTexture, nullptr, &textRectb);

		ss.str("");
		ss.clear();
		ss  << "h = "  << std::fixed << std::setprecision(2) << h;
		textSurface = TTF_RenderText_Solid(font, ss.str().c_str(), textColor);
		textTexture = SDL_CreateTextureFromSurface(renderer, textSurface);

		SDL_FreeSurface(textSurface);
		SDL_Rect textRecth = {0, 60, textSurface->w, textSurface->h};
		SDL_RenderCopy(renderer, textTexture, nullptr, &textRecth);
	}

	void drawVertices(){
        SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
        SDL_RenderClear(renderer);
        
        for(int i = 0; i != coordinates.size()-1; i++){
			SDL_SetRenderDrawColor(renderer, 150, 0, 205, 200);
			drawPoint(1, coordinates[i].x, coordinates[i].y);
        }
		drawText();
		
        SDL_RenderPresent(renderer);
    }

	void changeAttractor(int k = 100){
		switch(current_method){
			case Euler:
				setCoordinates(&Attractor::EulerMethod);
				break;
			case MiddlePoint:
				setCoordinates(&Attractor::middlePointMethod);
				break;
			case RK_4:
				setCoordinates(&Attractor::RK4);
				break;
			case Trapezoid:
				setCoordinates(&Attractor::reverseTrapezoidMethod);
				break;			
			case PredictorCorrector:
				setCoordinates(&Attractor::predictorCorrectorMethod);
				break;
			case DOPRI_8:
				setCoordinates(&Attractor::DOPRI8);
				break;
			default:
				break;
		}
		transformAttractor(k, true);
		for (auto& rotate: rotates){
			rotateAttractor(rotate);
		}
	}

	void drawAttractor() {
		int k = 100;
		int current_rotate = 10;
        SDL_Event event;
        bool running = true;
        transformAttractor(k, true);
        drawVertices();
        while (running) {
            while (SDL_PollEvent(&event)) {
                if (event.type == SDL_QUIT) running = false;
                if (event.type == SDL_KEYDOWN) {
                    switch (event.key.keysym.sym) {
                        case SDLK_UP:		//rotate up
							rotates.push_back(Rotate{X, current_rotate});
                            rotateAttractor(rotates.back());
                            drawVertices();
                            break;
                        case SDLK_DOWN:		//rotate down
                            rotates.push_back(Rotate{X, -current_rotate});
							rotateAttractor(rotates.back());							
                            drawVertices();
                            break;
                        case SDLK_LEFT:		//rotate left
                            rotates.push_back(Rotate{Y, -current_rotate});
							rotateAttractor(rotates.back());							
                            drawVertices();
                            break;
                        case SDLK_RIGHT:	//rotate right
							rotates.push_back(Rotate{Y, current_rotate});
                            rotateAttractor(rotates.back());
                            drawVertices();
                            break;
                        case SDLK_PAGEUP: 	//zoom
							k += 5;
							transformAttractor(k, false, k - 5);
							drawVertices();
							break;
						case SDLK_PAGEDOWN: //move away
							k -= 5;
							transformAttractor(k, false, k + 5);
							drawVertices();
							break;
						case SDLK_COMMA: 	// change h
							h -= 0.03;
							changeAttractor(k);
							drawVertices();
							break;
						case SDLK_PERIOD: 	// change h
							h += 0.03;
							changeAttractor(k);
							drawVertices();
							break;
						case SDLK_RIGHTBRACKET:  //change b
							b += 0.01;
							changeAttractor(k);
							drawVertices();
							break;
						case SDLK_LEFTBRACKET: //change b
							b -= 0.01;
							changeAttractor(k);
							drawVertices();
							break;
						case SDLK_e:	//change method to EulerMethod
							b = 0.19;
							h = 0.1;
							current_method = Euler;
							count_point = 3000000;
							changeAttractor(k);
							drawVertices();
							break;
						case SDLK_m:	//change method to MiddlePointMethod
							b = 0.19;
							h = 0.6;
							current_method = MiddlePoint;
							count_point = 3000000;
							changeAttractor(k);
							drawVertices();
							break;
						case SDLK_r:	//change method to RK4
							b = 0.19;
							h = 1.1;
							current_method = RK_4;
							count_point = 3000000;
							changeAttractor(k);
							drawVertices();
							break;
						case SDLK_t:	//change method to reverseTrapezoidMethod
							b = 0.18;
							h = -0.55;
							current_method = Trapezoid;
							count_point = 3000000;
							changeAttractor(k);
							drawVertices();
							break;
						case SDLK_p:	//change method to predictorCorrectorMethod
							b = 0.18;
							h = -0.78;
							current_method = PredictorCorrector;
							count_point = 3000000;
							changeAttractor(k);
							drawVertices();
							break;
						case SDLK_d:	//change method to DOPRI8
							b = 0.18;
							h = 0.7;
							current_method = DOPRI_8;
							count_point = 22000;
							changeAttractor(k);
							drawVertices();
						default:
                            break;
                    }
                }
            }
            

            SDL_Delay(10);
        }
    }

	~Attractor() {
        SDL_DestroyRenderer(renderer);
        SDL_DestroyWindow(window);
        SDL_Quit();
    }

};

double directedHausdorffDistance(vector<Vector>& A, vector<Vector>& B){
	double max_distance = 0;
	for (auto& a: A){
		double min_distance = numeric_limits<double>::max();
		for(auto& b: B){
			double current_distance = euclideanDistance(a, b);
			min_distance = min(min_distance, current_distance);
		} 
		max_distance = max(max_distance, min_distance);
	}
	return max_distance;
}

double HausdorffDistance(std::vector<Vector>& standart, vector<Vector>& new_way){
	double h_A_B = directedHausdorffDistance(standart, new_way);
	double h_B_A = directedHausdorffDistance(new_way, standart);

	return max(h_A_B, h_B_A);
}



int main(){
	float h = -0.78;
	Attractor atr (3000000, 0.18, h);
	atr.drawAttractor();
}