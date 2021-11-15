#include <iostream>
#include <vector>
#include <list>
#include <cmath>

using namespace std;


class vect{
public:
    double x;
    double y ;
    double z;
    vect(){
        x = 0.0;
        y = 0.0;
        z = 0.0;
    }
    vect(double x_, double y_, double z_){
        x = x_;
        y = y_;
        z = z_;
    }
    double get_length(){
        return sqrt(x * x + y * y + z * z);
    }
};


vect operator-(vect v1, vect v2){
    vect ans(v2.x - v1.x, v2.y - v1.y, v2.z - v1.z);
    return ans;
}
vect operator*(double a, vect v1){
    vect ans(a * v1.x, a * v1.y, a * v1.z);
    return ans;
}
vect operator*(vect v1, double a){
    vect ans(a * v1.x, a * v1.y, a * v1.z);
    return ans;
}

vect vectorProduct(vect v1, vect v2){
    vect ans;
    ans.x = v1.x * v2.y - v1.y * v2.x;
    ans.y = v1.y* v2.z - v1.z * v2.y;
    ans.z = v1.z * v2.x - v1.x * v2.z;
    return ans;
}

double dist(vect v1, vect v2){
    return (v2 - v1).get_length();
}

class triangle{
public:
    vect point1, point2, point3;
    triangle(){
        vect p0(0, 0, 0);
        point1 = p0;
        point2 = p0;
        point3 = p0;
    }
    triangle(vect point1_, vect point2_, vect point3_){
        point1 = point1_;
        point2 = point2_;
        point3 = point3_;
    }
    double get_area(){
        return abs(vectorProduct(point2 - point1, point3 - point1).get_length()) * 0.5;
    }

};
int main(){
//    vect p1(0, 0, 0);
//    vect p2(0, 0, 0);
//    vect p3(0, 0, 0);
//    triangle tr(p1, p2, p3);
    int N = 100;
    vector<vect> contour(N);
    for(int i = 0; i < N; i++){
        vect p(cos(i  / double(N) * 2 * 3.14159), sin(i  / double(N) * 2 * 3.14159), 0);
        contour[i] = p;
    }
    double mean_x = 0.0, mean_y = 0.0, mean_z = 0.0;
    for(int i=0; i < contour.size(); i++){
        mean_x += contour[i].x;
        mean_y += contour[i].y;
        mean_z += contour[i].z;
    }
    mean_x /= N;
    mean_y /= N;
    mean_z /= N;
    vect mean_v(mean_x, mean_y, mean_z);
    for(int i=0; i < N; i++){
        contour[i] = contour[i] - mean_v;
    }

    int TriagNMAX = 1e5;
    int pointsNMAX = 1e5;
    vector<triangle> triangles(TriagNMAX);
    vector<vect> points(pointsNMAX);
    vector<list<int> > pointTriangles(pointsNMAX);
    double Dl = 0.0;
    double meanDist = 0.0;
    double minDist = 1e10;
    for(int i=0; i<N; i++){
        Dl += dist(contour[i], contour[(i + 1) % N]);
        meanDist += contour[i].get_length();
        minDist = min(minDist, contour[i].get_length());
    }
    Dl /= N;
    meanDist /= N;
    cout << Dl << " " << meanDist << '\n';
    // first step - put contour into points array
    for(int i=0; i < N; i++){
        points[i] = contour[i];
    }
    // next steps
    int point_ind = 0, triangle_ind = 0;
    int cur_N = N;
    double alpha = 1.0;
    while(1 - Dl / meanDist > 0){
        cout << alpha << '\n';
        alpha = 1 - Dl / meanDist;
        int checkN = int(meanDist / Dl); // if i % checkN == 0 we need to delete this point
        int stop = point_ind + cur_N;
        int deleted = 0;
        int i0 = point_ind;
        cout << cur_N << '\n';
        while(point_ind < stop){
            if(point_ind % checkN != 0 || point_ind >= stop - 2){
                int point_ind2 = point_ind + 1;
                if(point_ind + 1 == stop)
                    point_ind2 = i0;
                points[point_ind + cur_N] = alpha * points[point_ind];
                points[point_ind2 + cur_N] = alpha * points[point_ind2];
                triangle tr1(points[point_ind], points[point_ind2], points[point_ind2 + cur_N]);
                triangle tr2(points[point_ind], points[point_ind + cur_N], points[point_ind2 + cur_N]);
                triangles[triangle_ind] = tr1;
                triangles[triangle_ind + 1] = tr2;
                pointTriangles[point_ind].push_back(triangle_ind);
                pointTriangles[point_ind].push_back(triangle_ind + 1);
                pointTriangles[point_ind2].push_back(triangle_ind + 1);
                pointTriangles[point_ind + cur_N].push_back(triangle_ind);
                pointTriangles[point_ind2 + cur_N].push_back(triangle_ind);
                pointTriangles[point_ind2 + cur_N].push_back(triangle_ind + 1);
                point_ind += 1;
                triangle_ind += 2;
            }
            else{
                deleted++;
                points[point_ind + cur_N] = alpha * points[point_ind];
                points[point_ind + cur_N + 1] = alpha * points[point_ind + 2];
                triangle tr1(points[point_ind], points[point_ind + 1], points[point_ind + cur_N]);
                triangle tr2(points[point_ind + cur_N], points[point_ind + 1], points[point_ind + cur_N + 1]);
                triangle tr3(points[point_ind], points[point_ind + 1], points[point_ind + cur_N]);
                pointTriangles[point_ind].push_back(triangle_ind);
                pointTriangles[point_ind + 1].push_back(triangle_ind);
                pointTriangles[point_ind + 1].push_back(triangle_ind + 1);
                pointTriangles[point_ind + 1].push_back(triangle_ind + 2);
                pointTriangles[point_ind + 2].push_back(triangle_ind + 2);
                pointTriangles[point_ind + cur_N].push_back(triangle_ind);
                pointTriangles[point_ind + cur_N].push_back(triangle_ind + 1);
                pointTriangles[point_ind + cur_N + 1].push_back(triangle_ind + 1);
                pointTriangles[point_ind + cur_N + 1].push_back(triangle_ind + 2);
                point_ind += 2;
                triangle_ind += 3;
            }
        }
        cur_N -= deleted;
        meanDist *= alpha;
    }
    //last step - connect all remain points with centre
    points[point_ind + cur_N] = vect(0, 0, 0);
    int stop = point_ind + cur_N, i0 = point_ind;
    while(point_ind < stop - 1){
        triangle tr(points[point_ind], points[point_ind + 1], points[stop]);
        triangles[triangle_ind] = tr;
        pointTriangles[point_ind].push_back(triangle_ind);
        pointTriangles[point_ind + 1].push_back(triangle_ind);
        pointTriangles[stop].push_back(triangle_ind);
        triangle_ind += 1;
        point_ind += 1;
    }
    triangle tr(points[point_ind], points[i0], points[stop]);
    triangles[triangle_ind] = tr;
    pointTriangles[point_ind].push_back(triangle_ind);
    pointTriangles[point_ind + 1].push_back(triangle_ind);
    pointTriangles[stop].push_back(triangle_ind);
    cout << triangle_ind << '\n';

     //first triangulation now is done

    return 0;
}