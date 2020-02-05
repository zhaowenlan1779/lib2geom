#include <2geom/2geom.h>


int main(){
    Geom::Rect rect1(0,0,1,1);
    Geom::Rect rect2(0.5,0.5,1.5,1.5);

    return rect1.intersects(rect2)? 0 : 1;
}
