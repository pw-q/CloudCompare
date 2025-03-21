//##########################################################################
//#                                                                        #
//#                              CLOUDCOMPARE                              #
//#                                                                        #
//#  This program is free software; you can redistribute it and/or modify  #
//#  it under the terms of the GNU General Public License as published by  #
//#  the Free Software Foundation; version 2 or later of the License.      #
//#                                                                        #
//#  This program is distributed in the hope that it will be useful,       #
//#  but WITHOUT ANY WARRANTY; without even the implied warranty of        #
//#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the          #
//#  GNU General Public License for more details.                          #
//#                                                                        #
//#          COPYRIGHT: EDF R&D / TELECOM ParisTech (ENST-TSI)             #
//#                                                                        #
//##########################################################################

#ifndef CC_TERRAIN_MODEL_TOOL_HEADER
#define CC_TERRAIN_MODEL_TOOL_HEADER

// Local
#include "cc2.5DimEditor.h"

// Qt
#include <QDialog>

class ccGenericPointCloud;
class ccPointCloud;
class ccPolyline;

namespace Ui {
class TerrainModelDialog;
}

//! Volume calculation tool (dialog)
class ccTerrainModelTool : public QDialog, public cc2Point5DimEditor {
  Q_OBJECT

public:
  //! Default constructor
  ccTerrainModelTool(ccGenericPointCloud *cloud, ccPolyline *edge = nullptr,
                     QWidget *parent = nullptr);

  //! Destructor
  ~ccTerrainModelTool();

  virtual double getGridStep() const override;
  virtual unsigned char getProjectionDimension() const override;
  virtual ccRasterGrid::ProjectionType getTypeOfProjection() const override;

  static bool
  generateMode(ccRasterGrid &grid, ccGenericPointCloud *ground,
               ccPolyline *edge, const ccBBox &gridBox, unsigned char vertDim,
               double gridStep, unsigned gridWidth, unsigned gridHeight,
               ccRasterGrid::ProjectionType projectionType,
               ccRasterGrid::EmptyCellFillOption groundEmptyCellFillStrategy,
               int krigingKNN, QWidget *parentWidget = nullptr);

protected:
  // Inherited from cc2Point5DimEditor
  virtual void gridIsUpToDate(bool state) override;

  //! Accepts the dialog and save settings
  void saveSettingsAndAccept();

  //! Save persistent settings and 'accept' dialog
  void saveSettings();

  //! Load persistent settings
  void loadSettings();

  void slotGenerate();

  void generateAndDisplay();

protected:
  Ui::TerrainModelDialog *m_ui;

  //! associated cloud
  ccGenericPointCloud *m_cloud;

  //包围盒
  ccPolyline *m_edge;
};

#endif // CC_VOLUME_CALC_TOOL_HEADER
