/* ParSGCT Code source file.
Copyright (c) 2015 Peter Edward Strazdins. All rights reserved.
Licensed under the terms of the BSD License as described in the LICENSE_SGCT file.
This comment must be retained in any redistributions of this source file.
*/

// 2D vector class supporting basic arithmetic
// written by Peter Strazdins, Jun 14

#ifndef VEC3D_INCLUDED
#define VEC3D_INCLUDED

#include <algorithm> // std::min
 
template <class T>
class Vec3D {
 public:
  union {
    struct {T x, y, z;};
    T v[3];
  };
  Vec3D<T> (T xv, T yv, T zv) {
    x = xv; y = yv; z = zv;
  }
  Vec3D<T> (T v) {
    x = y = z = v;
  }
  Vec3D<T> () {
    x = y = z = 0;
  }
 
  Vec3D<T> min(Vec3D<T> v) {
    return Vec3D<T> (std::min(x, v.x), std::min(y, v.y),  std::min(z, v.z));
  } 
  Vec3D<T> max(Vec3D<T> v) {
    return Vec3D<T> (std::max(x, v.x), std::max(y, v.y),  std::max(z, v.z));
  } 
  Vec3D<T> whereEq (Vec3D<T> v) {
    return Vec3D<T> (x == v.x, y == v.y, z == v.z);
  }

  T prod() {
    return x * y * z;
  }
  T sum() {
    return x + y + z;
  }

  Vec3D<T> operator + (const Vec3D<T>& v) const {  
    return Vec3D<T>  (x + v.x, y + v.y, z + v.z);
  }
  Vec3D<T> operator + (int v) const {  
    return Vec3D<T>  (x + v, y + v, z + v);
  }
  Vec3D<T> operator - (const Vec3D<T>& v) const {  
    return Vec3D<T>  (x - v.x, y - v.y, z - v.z);
  }
  Vec3D<T> operator - (int v) const {  
    return Vec3D<T>  (x - v, y - v, z - v);
  }
  Vec3D<T> operator * (const Vec3D<T>& v) const {  
    return Vec3D<T>  (x * v.x, y * v.y, z * v.z);
  }
  Vec3D<T> operator * (int v) const {  
    return Vec3D<T>  (x * v, y * v, z * v);
  }
  Vec3D<T> operator / (const Vec3D<T>& v) const {  
    return Vec3D<T>  (x / v.x, y / v.y, z / v.z);
  }
  friend Vec3D<T> operator / (int a, const Vec3D<T>& v) {  
    return Vec3D<T>  (a / v.x, a / v.y, a / v.z);
  }
  Vec3D<T> operator % (const Vec3D<T>& v) const {  
    return Vec3D<T>  (x % v.x, y % v.y, z % v.z);
  }
  friend bool operator == (const Vec3D<T>& u, const Vec3D<T>& v) {
    return (u.x == v.x && u.y == v.y && u.z == v.z);
  }
  friend bool operator <= (const Vec3D<T>& u, const Vec3D<T>& v) {
    return (u.x <= v.x && u.y <= v.y && u.z <= v.z);
  }  
  
};

#define V3DFMT "%d,%d,%d"
#define V3DFMTSZ "%dx%dx%d"
#define V3DLST(v) (v).x,(v).y,(v).z

#endif /*VEC3D_INCLUDED*/
