//##########################################################################
//#                                                                        #
//#                    ZOOMLION PLUGIN: ccCompass                      #
//#                                                                        #
//#  This program is free software; you can redistribute it and/or modify  #
//#  it under the terms of the GNU General Public License as published by  #
//#  the Free Software Foundation; version 2 of the License.               #
//#                                                                        #
//#  This program is distributed in the hope that it will be useful,       #
//#  but WITHOUT ANY WARRANTY; without even the implied warranty of        #
//#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         #
//#  GNU General Public License for more details.                          #
//#                                                                        #
//#                     COPYRIGHT: Sam Thiele  2017                        #
//#                                                                        #
//##########################################################################

#include "ccLineation.h"

// pass ctors straight to PointPair
ccLineation::ccLineation(ccPointCloud *associatedCloud)
    : ccPointPair(associatedCloud) {
  updateMetadata();
}

ccLineation::ccLineation(ccPolyline *obj) : ccPointPair(obj) {
  updateMetadata();
}

void ccLineation::updateMetadata() {
  QVariantMap map;

  // add metadata tag defining the ccCompass class type
  map.insert("ccCompassType", "Lineation");

  // calculate trace orientation (trend/plunge)
  if (size() ==
      2) // can't calculate orientation of something smaller than this...
  {
    CCVector3f dir = getDirection();
    dir.normalize();
    float trend = 0.0f;
    float plunge = 0.0f;

    if (dir.x + dir.y + dir.z == 0) // special case: null direction
    {
      trend = 0;
      plunge = 0;
    } else if (dir.z > 0.9999999 ||
               dir.z < -0.9999999) // special case: dir is vertical (= 0,0,1 or
                                   // 0,0,-1)
    {
      trend = 0;
      if (dir.z < 0)
        plunge = 90;
      else
        plunge = -90;
    } else // normal cases...
    {
      CCVector3f hzComp = CCVector3f(dir.x, dir.y, 0);
      hzComp.normalize();

      // calculate plunge: plunge = angle between vector & vector projected onto
      // horizontal (x,y) plane
      plunge = std::acos(dir.dot(hzComp)) *
               (180 / M_PI); // plunge measured from horizontal (in degrees)
      if (dir.z > 0) // lineations pointing towards the sky have negative
                     // plunges
        plunge *= -1;

      // calculate trend (N.B. I have very little idea how exactly this code
      // work, it's kinda magic) [c.f.
      //http://stackoverflow.com/questions/14066933/direct-way-of-computing-clockwise-angle-between-2-vectors
      //]
      CCVector3f N(0, 1, 0); // north vector
      float dot = hzComp.dot(N);
      float det = CCVector3f(0, 0, 1).dot(hzComp.cross(N));
      trend =
          std::atan2(det, dot) *
          (180 / M_PI); // heading measured clockwise from north (in degrees)
      if (trend < 0)
        trend += 360;
    }

    // store trend and plunge info
    // CCVector3d Pg = cloud->toGlobal3d(*P);
    CCVector3d s = toGlobal3d(*getPoint(0)); // start point
    CCVector3d e = toGlobal3d(*getPoint(1)); // end point
    float length = (s - e).norm();

    map.insert("Sx", s.x);
    map.insert("Sy", s.y);
    map.insert("Sz", s.z);
    map.insert("Ex", e.x);
    map.insert("Ey", e.y);
    map.insert("Ez", e.z);
    map.insert("Trend", trend);
    map.insert("Plunge", plunge);
    map.insert("Length", length * getGlobalScale());

    // store metadata
    setMetaData(map, true);

    // update name
    QString lengthstr = QString("").asprintf("%.1f on ", length);
    QString trendAndPlungeStr = QString("%2->%3")
                                    .arg((int)plunge, 2, 10, QChar('0'))
                                    .arg((int)trend, 3, 10, QChar('0'));
    QString namestr = lengthstr + trendAndPlungeStr;
    setName(namestr);
  }
}

// returns true if object is a lineation
bool ccLineation::isLineation(ccHObject *object) {
  if (object->hasMetaData("ccCompassType")) {
    return object->getMetaData("ccCompassType")
        .toString()
        .contains("Lineation");
  }
  return false;
}
