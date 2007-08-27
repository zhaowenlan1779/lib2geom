#ifndef __2GEOM_SHAPE_H
#define __2GEOM_SHAPE_H

#include <vector>
#include <set>

#include "region.h"

//TODO: BBOX optimizations

namespace Geom {

enum {
  BOOLOP_JUST_A  = 1,
  BOOLOP_JUST_B  = 2,
  BOOLOP_BOTH    = 4,
  BOOLOP_NEITHER = 8
};

enum {
  BOOLOP_NULL         = 0,
  BOOLOP_INTERSECT    = BOOLOP_BOTH,
  BOOLOP_SUBTRACT_A_B = BOOLOP_JUST_B,
  BOOLOP_IDENTITY_A   = BOOLOP_JUST_A | BOOLOP_BOTH,
  BOOLOP_SUBTRACT_B_A = BOOLOP_JUST_A,
  BOOLOP_IDENTITY_B   = BOOLOP_JUST_B | BOOLOP_BOTH,
  BOOLOP_EXCLUSION    = BOOLOP_JUST_A | BOOLOP_JUST_B,
  BOOLOP_UNION        = BOOLOP_JUST_A | BOOLOP_JUST_B | BOOLOP_BOTH
};

class Shape {
    Regions content;
    mutable bool fill;
    //friend Shape shape_region_boolean(bool rev, Shape const & a, Region const & b);
    friend CrossingSet crossings_between(Shape const &a, Shape const &b);
    friend Shape shape_boolean(bool rev, Shape const &, Shape const &, CrossingSet const &);
    friend Shape shape_boolean(Shape const &a, Shape const &b, unsigned);
    friend Shape shape_boolean(Shape const &a, Shape const &b, unsigned, CrossingSet const &);
    std::vector<Rect> bounds() const;

  public:
    Shape() : fill(true) {}
    explicit Shape(Region const & r) {
        content = Regions(1, r);
        fill = r.fill;
    }
    explicit Shape(Regions const & r) : content(r) { update_fill(); }
    explicit Shape(bool f) : fill(f) {}
    Shape(Regions const & r, bool f) : content(r), fill(f) {}
    
    Regions getContent() const { return content; }
    bool isFill() const { return fill; }
    
    unsigned size() const { return content.size(); }
    const Region &operator[](unsigned ix) const { return content[ix]; }
    
    Shape inverse() const;
    Shape operator*(Matrix const &m) const;
    
    bool contains(Point const &p) const;
    
    bool inside_invariants() const;  //semi-slow & easy to violate : checks that the insides are inside, the outsides are outside
    bool region_invariants() const; //semi-slow                    : checks for self crossing
    bool cross_invariants() const; //slow                          : checks that everything is disjoint
    bool invariants() const;      //vera slow (combo rombo, checks the above)

  private:     
    void update_fill() const {
        unsigned ix = outer_index(content);
        if(ix < size())
            fill = content[ix].fill;
        else if(size() > 0)
            fill = content.front().fill;
        else
            fill = true;
    }
};

CrossingSet crossings_between(Shape const &a, Shape const &b);

Shape shape_boolean(bool rev, Shape const &, Shape const &, CrossingSet const &);
Shape shape_boolean(bool rev, Shape const &, Shape const &);

Shape shape_boolean(Shape const &, Shape const &, unsigned flags);
Shape shape_boolean(Shape const &, Shape const &, unsigned flags, CrossingSet &);

Shape sanitize_paths(std::vector<Path> ps);

}

#endif
