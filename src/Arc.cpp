#include "Arc.h"

using namespace std;

ostream& operator<<(ostream& os, const Arc& a) {
    os<<"["<<a.i<<","<<a.j<<"]";
    return os;
}

