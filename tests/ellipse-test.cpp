/** @file
 * @brief Unit tests for Ellipse and related functions
 * Uses the Google Testing Framework
 *//*
 * Authors:
 *   Krzysztof Kosiński <tweenk.pl@gmail.com>
 *
 * Copyright 2015 Authors
 *
 * This library is free software; you can redistribute it and/or
 * modify it either under the terms of the GNU Lesser General Public
 * License version 2.1 as published by the Free Software Foundation
 * (the "LGPL") or, at your option, under the terms of the Mozilla
 * Public License Version 1.1 (the "MPL"). If you do not alter this
 * notice, a recipient may use your version of this file under either
 * the MPL or the LGPL.
 *
 * You should have received a copy of the LGPL along with this library
 * in the file COPYING-LGPL-2.1; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 * You should have received a copy of the MPL along with this library
 * in the file COPYING-MPL-1.1
 *
 * The contents of this file are subject to the Mozilla Public License
 * Version 1.1 (the "License"); you may not use this file except in
 * compliance with the License. You may obtain a copy of the License at
 * http://www.mozilla.org/MPL/
 *
 * This software is distributed on an "AS IS" basis, WITHOUT WARRANTY
 * OF ANY KIND, either express or implied. See the LGPL or the MPL for
 * the specific language governing rights and limitations.
 */

#include <iostream>
#include <glib.h>

#include <2geom/angle.h>
#include <2geom/ellipse.h>
#include <2geom/elliptical-arc.h>
#include <memory>

#include "testing.h"

#ifndef M_SQRT2
#  define M_SQRT2 1.41421356237309504880
#endif

using namespace Geom;

TEST(EllipseTest, Arcs) {
    Ellipse e(Point(5,10), Point(5, 10), 0);

    std::unique_ptr<EllipticalArc> arc1(e.arc(Point(5,0), Point(0,0), Point(0,10)));

    EXPECT_EQ(arc1->initialPoint(), Point(5,0));
    EXPECT_EQ(arc1->finalPoint(), Point(0,10));
    EXPECT_EQ(arc1->boundsExact(), Rect::from_xywh(0,0,5,10));
    EXPECT_EQ(arc1->center(), e.center());
    EXPECT_EQ(arc1->largeArc(), false);
    EXPECT_EQ(arc1->sweep(), false);

    std::unique_ptr<EllipticalArc> arc1r(e.arc(Point(0,10), Point(0,0), Point(5,0)));

    EXPECT_EQ(arc1r->boundsExact(), arc1->boundsExact());
    EXPECT_EQ(arc1r->sweep(), true);
    EXPECT_EQ(arc1r->largeArc(), false);

    std::unique_ptr<EllipticalArc> arc2(e.arc(Point(5,0), Point(10,20), Point(0,10)));

    EXPECT_EQ(arc2->boundsExact(), Rect::from_xywh(0,0,10,20));
    EXPECT_EQ(arc2->largeArc(), true);
    EXPECT_EQ(arc2->sweep(), true);

    std::unique_ptr<EllipticalArc> arc2r(e.arc(Point(0,10), Point(10,20), Point(5,0)));

    EXPECT_EQ(arc2r->boundsExact(), arc2->boundsExact());
    EXPECT_EQ(arc2r->largeArc(), true);
    EXPECT_EQ(arc2r->sweep(), false);

    // exactly half arc
    std::unique_ptr<EllipticalArc> arc3(e.arc(Point(5,0), Point(0,10), Point(5,20)));

    EXPECT_EQ(arc3->boundsExact(), Rect::from_xywh(0,0,5,20));
    EXPECT_EQ(arc3->largeArc(), false);
    EXPECT_EQ(arc3->sweep(), false);

    // inner point exactly at midpoint between endpoints
    std::unique_ptr<EllipticalArc> arc4(e.arc(Point(5,0), Point(2.5,5), Point(0,10)));

    EXPECT_EQ(arc4->initialPoint(), Point(5,0));
    EXPECT_EQ(arc4->finalPoint(), Point(0,10));
    EXPECT_EQ(arc4->boundsExact(), Rect::from_xywh(0,0,5,10));
    EXPECT_EQ(arc4->largeArc(), false);
    EXPECT_EQ(arc4->sweep(), false);

    std::unique_ptr<EllipticalArc> arc4r(e.arc(Point(0,10), Point(2.5,5), Point(5,0)));

    EXPECT_EQ(arc4r->initialPoint(), Point(0,10));
    EXPECT_EQ(arc4r->finalPoint(), Point(5,0));
    EXPECT_EQ(arc4r->boundsExact(), Rect::from_xywh(0,0,5,10));
    EXPECT_EQ(arc4r->largeArc(), false);
    EXPECT_EQ(arc4r->sweep(), true);
}

TEST(EllipseTest, AreNear) {
    Ellipse e1(Point(5.000001,10), Point(5,10), Angle::from_degrees(45));
    Ellipse e2(Point(5.000000,10), Point(5,10), Angle::from_degrees(225));
    Ellipse e3(Point(4.999999,10), Point(10,5), Angle::from_degrees(135));
    Ellipse e4(Point(5.000001,10), Point(10,5), Angle::from_degrees(315));

    EXPECT_TRUE(are_near(e1, e2, 1e-5));
    EXPECT_TRUE(are_near(e1, e3, 1e-5));
    EXPECT_TRUE(are_near(e1, e4, 1e-5));

    Ellipse c1(Point(20.000001,35.000001), Point(5.000001,4.999999), Angle::from_degrees(180.00001));
    Ellipse c2(Point(19.999999,34.999999), Point(4.999999,5.000001), Angle::from_degrees(179.99999));
    //std::cout << c1 << "\n" << c2 << std::endl;
    EXPECT_TRUE(are_near(c1, c2, 2e-5));

    EXPECT_FALSE(are_near(c1, e1, 1e-5));
    EXPECT_FALSE(are_near(c2, e1, 1e-5));
    EXPECT_FALSE(are_near(c1, e2, 1e-5));
    EXPECT_FALSE(are_near(c2, e2, 1e-5));
    EXPECT_FALSE(are_near(c1, e3, 1e-5));
    EXPECT_FALSE(are_near(c2, e3, 1e-5));
    EXPECT_FALSE(are_near(c1, e4, 1e-5));
    EXPECT_FALSE(are_near(c2, e4, 1e-5));
}

TEST(EllipseTest, Transformations) {
    Ellipse e(Point(5,10), Point(5,10), Angle::from_degrees(45));

    Ellipse er = e * Rotate::around(Point(5,10), Angle::from_degrees(45));
    Ellipse ercmp(Point(5,10), Point(5,10), Angle::from_degrees(90));
    //std::cout << e << "\n" << er << "\n" << ercmp << std::endl;
    EXPECT_TRUE(are_near(er, ercmp, 1e-12));

    Ellipse eflip = e * Affine(Scale(-1,1));
    Ellipse eflipcmp(Point(-5, 10), Point(5,10), Angle::from_degrees(135));
    EXPECT_TRUE(are_near(eflip, eflipcmp, 1e-12));
}

TEST(EllipseTest, TimeAt) {
    Ellipse e(Point(4, 17), Point(22, 34), 2);

    for (unsigned i = 0; i < 100; ++i) {
        Coord t = g_random_double_range(0, 2*M_PI);
        Point p = e.pointAt(t);
        Coord t2 = e.timeAt(p);
        EXPECT_FLOAT_EQ(t, t2);
    }
}

TEST(EllipseTest, LineIntersection) {
    Ellipse e(Point(0, 0), Point(3, 2), 0);
    Line l(Point(0, -2), Point(1, 0));

    std::vector<ShapeIntersection> xs = e.intersect(l);

    ASSERT_EQ(xs.size(), 2ul);

    // due to numeric imprecision when evaluating Ellipse,
    // the points may deviate by around 2e-16
    EXPECT_NEAR(xs[0].point()[X], 0, 1e-15);
    EXPECT_NEAR(xs[0].point()[Y], -2, 1e-15);
    EXPECT_NEAR(xs[1].point()[X], 9./5, 1e-15);
    EXPECT_NEAR(xs[1].point()[Y], 8./5, 1e-15);

    EXPECT_intersections_valid(e, l, xs, 1e-15);

    // Test with a degenerate ellipse
    auto degen = Ellipse({0, 0}, {3, 2}, 0);
    degen *= Scale(1.0, 0.0); // Squash to the X-axis interval [-3, 3].

    g_random_set_seed(0xCAFECAFE);
    // Intersect with a line
    for (size_t _ = 0; _ < 10'000; _++) {
        auto line = Line(Point(g_random_double_range(-3.0, 3.0), g_random_double_range(-3.0, -1.0)),
                         Point(g_random_double_range(-3.0, 3.0), g_random_double_range(1.0, 3.0)));
        auto xings = degen.intersect(line);
        EXPECT_EQ(xings.size(), 2u);
        EXPECT_intersections_valid(degen, line, xings, 1e-14);
    }
    // Intersect with another, non-degenerate ellipse
    for (size_t _ = 0; _ < 10'000; _++) {
        auto other = Ellipse(Point(g_random_double_range(-1.0, 1.0), g_random_double_range(-1.0, 1.0)),
                             Point(g_random_double_range(1.0, 2.0), g_random_double_range(1.0, 3.0)), 0);
        auto xings = degen.intersect(other);
        EXPECT_intersections_valid(degen, other, xings, 1e-14);
    }
    // Intersect with another ellipse which is also degenerate
    for (size_t _ = 0; _ < 10'000; _++) {
        auto other = Ellipse({0, 0}, {1, 1}, 0); // Unit circle
        other *= Scale(0.0, g_random_double_range(0.5, 4.0)); // Squash to Y axis
        other *= Rotate(g_random_double_range(-1.5, 1.5)); // Rotate a little (still passes through the origin)
        other *= Translate(g_random_double_range(-2.9, 2.9), 0.0);
        auto xings = degen.intersect(other);
        EXPECT_EQ(xings.size(), 4u);
        EXPECT_intersections_valid(degen, other, xings, 1e-14);
    }
}

TEST(EllipseTest, EllipseIntersection) {
    Ellipse e1;
    Ellipse e2;
    std::vector<ShapeIntersection> xs;

    e1.set(Point(300, 300), Point(212, 70), -0.785);
    e2.set(Point(250, 300), Point(230, 90), 1.321);
    xs = e1.intersect(e2);
    EXPECT_EQ(xs.size(), 4ul);
    EXPECT_intersections_valid(e1, e2, xs, 4e-10);

    e1.set(Point(0, 0), Point(1, 1), 0);
    e2.set(Point(0, 1), Point(1, 1), 0);
    xs = e1.intersect(e2);
    EXPECT_EQ(xs.size(), 2ul);
    EXPECT_intersections_valid(e1, e2, xs, 1e-10);

    e1.set(Point(0, 0), Point(1, 1), 0);
    e2.set(Point(1, 0), Point(1, 1), 0);
    xs = e1.intersect(e2);
    EXPECT_EQ(xs.size(), 2ul);
    EXPECT_intersections_valid(e1, e2, xs, 1e-10);

    // === Test detection of external tangency between ellipses ===
    // Perpendicular major axes
    e1.set({0, 0}, {5, 3}, 0); // rightmost point (5, 0)
    e2.set({6, 0}, {1, 2}, 0); // leftmost point (5, 0)
    xs = e1.intersect(e2);
    ASSERT_GT(xs.size(), 0);
    EXPECT_intersections_valid(e1, e2, xs, 1e-10);
    EXPECT_TRUE(are_near(xs[0].point(), Point(5, 0)));

    // Collinear major axes
    e1.set({30, 0}, {9, 1}, 0); // leftmost point (21, 0)
    e2.set({18, 0}, {3, 2}, 0); // rightmost point (21, 0)
    xs = e1.intersect(e2);
    ASSERT_GT(xs.size(), 0);
    EXPECT_intersections_valid(e1, e2, xs, 1e-10);
    EXPECT_TRUE(are_near(xs[0].point(), Point(21, 0)));

    // Circles not aligned to an axis (Pythagorean triple: 3^2 + 4^2 == 5^2)
    e1.set({0, 0}, {3, 3}, 0); // radius 3
    e2.set({3, 4}, {2, 2}, 0); // radius 2
    // We know 2 + 3 == 5 == distance((0, 0), (3, 4)) so there's an external tangency
    // between these circles, at a point at distance 3 from the origin, on the line x = 0.75 y.
    xs = e1.intersect(e2);
    ASSERT_GT(xs.size(), 0);
    EXPECT_intersections_valid(e1, e2, xs, 1e-6);

    // === Test the detection of internal tangency between ellipses ===
    // Perpendicular major axes
    e1.set({0, 0}, {8, 17}, 0); // rightmost point (8, 0)
    e2.set({6, 0}, {2, 1}, 0); // rightmost point (8, 0)
    xs = e1.intersect(e2);
    ASSERT_GT(xs.size(), 0);
    EXPECT_intersections_valid(e1, e2, xs, 1e-10);
    EXPECT_TRUE(are_near(xs[0].point(), Point(8, 0)));

    // Collinear major axes
    e1.set({30, 0}, {9, 5}, 0); // rightmost point (39, 0)
    e2.set({36, 0}, {3, 1}, 0); // rightmost point (39, 0)
    xs = e1.intersect(e2);
    ASSERT_GT(xs.size(), 0);
    EXPECT_intersections_valid(e1, e2, xs, 1e-6);
    EXPECT_TRUE(are_near(xs[0].point(), Point(39, 0)));

    // Circles not aligned to an axis (Pythagorean triple: 3^2 + 4^2 == 5^2)
    e1.set({4, 3}, {5, 5}, 0); // Passes through (0, 0), center on the line y = 0.75 x
    e2.set({8, 6}, {10, 10}, 0); // Also passes through (0, 0), center on the same line.
    xs = e1.intersect(e2);
    ASSERT_GT(xs.size(), 0);
    EXPECT_intersections_valid(e1, e2, xs, 1e-6);
    EXPECT_TRUE(are_near(xs[0].point(), Point(0, 0)));
}

TEST(EllipseTest, BezierIntersection) {
    Ellipse e(Point(300, 300), Point(212, 70), -3.926);
    D2<Bezier> b(Bezier(100, 300, 100, 500), Bezier(100, 100, 500, 500));

    std::vector<ShapeIntersection> xs = e.intersect(b);

    EXPECT_EQ(xs.size(), 2ul);
    EXPECT_intersections_valid(e, b, xs, 6e-12);
}

TEST(EllipseTest, Coefficients) {
    std::vector<Ellipse> es;
    es.emplace_back(Point(-15,25), Point(10,15), Angle::from_degrees(45).radians0());
    es.emplace_back(Point(-10,33), Point(40,20), M_PI);
    es.emplace_back(Point(10,-33), Point(40,20), Angle::from_degrees(135).radians0());
    es.emplace_back(Point(-10,-33), Point(50,10), Angle::from_degrees(330).radians0());

    for (auto & i : es) {
        Coord a, b, c, d, e, f;
        i.coefficients(a, b, c, d, e, f);
        Ellipse te(a, b, c, d, e, f);
        EXPECT_near(i, te, 1e-10);
        for (Coord t = -5; t < 5; t += 0.125) {
            Point p = i.pointAt(t);
            Coord eq = a*p[X]*p[X] + b*p[X]*p[Y] + c*p[Y]*p[Y]
              + d*p[X] + e*p[Y] + f;
            EXPECT_NEAR(eq, 0, 1e-10);
        }
    }
}

TEST(EllipseTest, UnitCircleTransform) {
    std::vector<Ellipse> es;
    es.emplace_back(Point(-15,25), Point(10,15), Angle::from_degrees(45));
    es.emplace_back(Point(-10,33), Point(40,20), M_PI);
    es.emplace_back(Point(10,-33), Point(40,20), Angle::from_degrees(135));
    es.emplace_back(Point(-10,-33), Point(50,10), Angle::from_degrees(330));

    for (auto & e : es) {
        EXPECT_near(e.unitCircleTransform() * e.inverseUnitCircleTransform(), Affine::identity(), 1e-8);

        for (Coord t = -1; t < 10; t += 0.25) {
            Point p = e.pointAt(t);
            p *= e.inverseUnitCircleTransform();
            EXPECT_near(p.length(), 1., 1e-10);
            p *= e.unitCircleTransform();
            EXPECT_near(e.pointAt(t), p, 1e-10);
        }
    }
}

TEST(EllipseTest, PointAt) {
    Ellipse a(Point(0,0), Point(10,20), 0);
    EXPECT_near(a.pointAt(0), Point(10,0), 1e-10);
    EXPECT_near(a.pointAt(M_PI/2), Point(0,20), 1e-10);
    EXPECT_near(a.pointAt(M_PI), Point(-10,0), 1e-10);
    EXPECT_near(a.pointAt(3*M_PI/2), Point(0,-20), 1e-10);

    Ellipse b(Point(0,0), Point(10,20), M_PI/2);
    EXPECT_near(b.pointAt(0), Point(0,10), 1e-10);
    EXPECT_near(b.pointAt(M_PI/2), Point(-20,0), 1e-10);
    EXPECT_near(b.pointAt(M_PI), Point(0,-10), 1e-10);
    EXPECT_near(b.pointAt(3*M_PI/2), Point(20,0), 1e-10);
}

TEST(EllipseTest, UnitTangentAt) {
    Ellipse a(Point(14,-7), Point(20,10), 0);
    Ellipse b(Point(-77,23), Point(40,10), Angle::from_degrees(45));

    EXPECT_near(a.unitTangentAt(0), Point(0,1), 1e-12);
    EXPECT_near(a.unitTangentAt(M_PI/2), Point(-1,0), 1e-12);
    EXPECT_near(a.unitTangentAt(M_PI), Point(0,-1), 1e-12);
    EXPECT_near(a.unitTangentAt(3*M_PI/2), Point(1,0), 1e-12);

    EXPECT_near(b.unitTangentAt(0), Point(-M_SQRT2/2, M_SQRT2/2), 1e-12);
    EXPECT_near(b.unitTangentAt(M_PI/2), Point(-M_SQRT2/2, -M_SQRT2/2), 1e-12);
    EXPECT_near(b.unitTangentAt(M_PI), Point(M_SQRT2/2, -M_SQRT2/2), 1e-12);
    EXPECT_near(b.unitTangentAt(3*M_PI/2), Point(M_SQRT2/2, M_SQRT2/2), 1e-12);
}

TEST(EllipseTest, Bounds)
{
    // Create example ellipses
    std::vector<Ellipse> es;
    es.emplace_back(Point(-15,25), Point(10,15), Angle::from_degrees(45));
    es.emplace_back(Point(-10,33), Point(40,20), M_PI);
    es.emplace_back(Point(10,-33), Point(40,20), Angle::from_degrees(111));
    es.emplace_back(Point(-10,-33), Point(50,10), Angle::from_degrees(222));

    // for reproducibility
    g_random_set_seed(1234);

    for (auto & e : es) {
        Rect r = e.boundsExact();
        Rect f = e.boundsFast();
        for (unsigned j = 0; j < 10000; ++j) {
            Coord t = g_random_double_range(-M_PI, M_PI);
            auto const p = e.pointAt(t);
            EXPECT_TRUE(r.contains(p));
            EXPECT_TRUE(f.contains(p));
        }
    }

    Ellipse e(Point(0,0), Point(10, 10), M_PI);
    Rect bounds = e.boundsExact();
    Rect coarse = e.boundsFast();
    EXPECT_EQ(bounds, Rect(Point(-10,-10), Point(10,10)));
    EXPECT_TRUE(bounds.contains(e.pointAt(0)));
    EXPECT_TRUE(bounds.contains(e.pointAt(M_PI/2)));
    EXPECT_TRUE(bounds.contains(e.pointAt(M_PI)));
    EXPECT_TRUE(bounds.contains(e.pointAt(3*M_PI/2)));
    EXPECT_TRUE(bounds.contains(e.pointAt(2*M_PI)));
    EXPECT_TRUE(coarse.contains(e.pointAt(0)));
    EXPECT_TRUE(coarse.contains(e.pointAt(M_PI/2)));
    EXPECT_TRUE(coarse.contains(e.pointAt(M_PI)));
    EXPECT_TRUE(coarse.contains(e.pointAt(3*M_PI/2)));
    EXPECT_TRUE(coarse.contains(e.pointAt(2*M_PI)));

    e = Ellipse(Point(0,0), Point(10, 10), M_PI/2);
    bounds = e.boundsExact();
    coarse = e.boundsFast();
    EXPECT_EQ(bounds, Rect(Point(-10,-10), Point(10,10)));
    EXPECT_TRUE(bounds.contains(e.pointAt(0)));
    EXPECT_TRUE(bounds.contains(e.pointAt(M_PI/2)));
    EXPECT_TRUE(bounds.contains(e.pointAt(M_PI)));
    EXPECT_TRUE(bounds.contains(e.pointAt(3*M_PI/2)));
    EXPECT_TRUE(bounds.contains(e.pointAt(2*M_PI)));
    EXPECT_TRUE(coarse.contains(e.pointAt(0)));
    EXPECT_TRUE(coarse.contains(e.pointAt(M_PI/2)));
    EXPECT_TRUE(coarse.contains(e.pointAt(M_PI)));
    EXPECT_TRUE(coarse.contains(e.pointAt(3*M_PI/2)));
    EXPECT_TRUE(coarse.contains(e.pointAt(2*M_PI)));
}
