#include <stdio.h>
#include <SDL2/SDL.h>
#include <iostream>
#include <cmath>
#include <algorithm> 
#include <vector>
#include <GL/gl.h>
#include <SDL2/SDL_ttf.h>
#include <string>
#include <iomanip>
#include <sstream>
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
};

class Attractor{
public:
	int count_point;
	double b;
	vector <Vector> coordinates;
	double h = 0.6;
	//double h = 0.5 // 0.9 // 0.6;
	enum Orient {X, Y, Z};

	struct Rotate{
		Orient orient;
		int angle;
	};

	vector <Rotate> rotates;

	int height = 1000;
    int width = 1000;
    SDL_Window* window;
	SDL_Texture* textTexture;
	SDL_Surface* textSurface;
    SDL_Renderer* renderer;
	TTF_Font* font;

	Attractor (int count_point, double b): count_point(count_point), b(b) {
		setСoordinates();

		SDL_Init(SDL_INIT_VIDEO);
        window = SDL_CreateWindow("Attractor",
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
		
		font = TTF_OpenFont("/home/sun/Документы/VblchMat/Roboto/Roboto-VariableFont_wdth,wght.ttf", 24); // 24 — размер шрифта
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

	void setСoordinates(){
		coordinates.clear();
		coordinates.push_back(Vector(0.1, 0.2, 0.2));
		for (int i = 0; i < count_point; i++){
			if (i%10 == 0){
				coordinates.push_back(RK4(coordinates.back()));
			}
		}
		for (int i = 0; i < 100; i++){
			cout << coordinates[i].x << ' ' << coordinates[i].y << ' ' << coordinates[i].z << '\n';
		}
	}

	Vector eulerMethod(Vector point){
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
		SDL_Color textColor = {255, 255, 255, 255}; // Белый цвет (R, G, B, A)

		stringstream ss;
		ss  << "b = " << std::fixed << std::setprecision(2) << b;
		textSurface = TTF_RenderText_Solid(font, ss.str().c_str(), textColor);
		textTexture = SDL_CreateTextureFromSurface(renderer, textSurface);

		SDL_FreeSurface(textSurface); // Освобождаем поверхность, так как она больше не нужна
		SDL_Rect textRectb = {0, 0, textSurface->w, textSurface->h}; // Позиция и размер текста
		SDL_RenderCopy(renderer, textTexture, nullptr, &textRectb);

		ss.str("");
		ss.clear();
		ss  << "h = "  << std::fixed << std::setprecision(2) << h;
		textSurface = TTF_RenderText_Solid(font, ss.str().c_str(), textColor);
		textTexture = SDL_CreateTextureFromSurface(renderer, textSurface);

		SDL_FreeSurface(textSurface); // Освобождаем поверхность, так как она больше не нужна
		SDL_Rect textRecth = {0, 30, textSurface->w, textSurface->h}; // Позиция и размер текста
		SDL_RenderCopy(renderer, textTexture, nullptr, &textRecth);
	}

	void drawVertices(){
        SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
        SDL_RenderClear(renderer);
        
        for(int i = 0; i != coordinates.size()-1; i++){
			SDL_SetRenderDrawColor(renderer, 150, 0, 205, 200);
			drawPoint(1, coordinates[i].x, coordinates[i].y);
        }
        //SDL_Delay(10);
		drawText();
		
        SDL_RenderPresent(renderer);
    }

	void changeAttractor(int k){
		setСoordinates();
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
                        case SDLK_UP:
							rotates.push_back(Rotate{X, current_rotate});
                            rotateAttractor(rotates.back());
                            drawVertices();
                            break;
                        case SDLK_DOWN:
                            rotates.push_back(Rotate{X, -current_rotate});
							rotateAttractor(rotates.back());							
                            drawVertices();
                            break;
                        case SDLK_LEFT:
                            rotates.push_back(Rotate{Y, -current_rotate});
							rotateAttractor(rotates.back());							
                            drawVertices();
                            break;
                        case SDLK_RIGHT:
							rotates.push_back(Rotate{Y, current_rotate});
                            rotateAttractor(rotates.back());
                            drawVertices();
                            break;
                        case SDLK_PAGEUP:
							k += 5;
							transformAttractor(k, false, k - 5);
							drawVertices();
							break;
						case SDLK_PAGEDOWN:
							k -= 5;
							transformAttractor(k, false, k + 5);
							drawVertices();
							break;
						case SDLK_COMMA: // change h
							h -= 0.1;
							changeAttractor(k);
							drawVertices();
							break;
						case SDLK_PERIOD: // change h
							h += 0.1;
							changeAttractor(k);
							drawVertices();
							break;
						case SDLK_RIGHTBRACKET:
							b += 0.01;
							changeAttractor(k);
							drawVertices();
							break;
						case SDLK_LEFTBRACKET:
							b -= 0.01;
							changeAttractor(k);
							drawVertices();
							break;
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

int main(){
	Attractor atr (3000000, 0.19);
	//Attractor atr (3000000, 0.19). h = 0.6;
	atr.drawAttractor();
}