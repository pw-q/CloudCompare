//##########################################################################
//#                                                                        #
//#                              ZOOMLION                              #
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

#include "ccTerrainModel.h"
#include "ccPersistentSettings.h"
#include "ui_terrainModelDlg.h"

// Local

// qCC_db
#include <ccPointCloud.h>
#include <ccPolyline.h>
#include <ccProgressDialog.h>
#include <ccScalarField.h>

// qCC_gl
#include <ccGLWindowInterface.h>

// Qt
#include <QClipboard>
#include <QMessageBox>
#include <QSettings>

// System
#include <cassert>

ccTerrainModelTool::ccTerrainModelTool(ccGenericPointCloud *cloud,
                                       ccPolyline *edge /*=nullptr*/,
                                       QWidget *parent /*=nullptr*/)
    : QDialog(parent, Qt::WindowMaximizeButtonHint | Qt::WindowCloseButtonHint),
      m_cloud(cloud), m_edge(edge), cc2Point5DimEditor(),
      m_ui(new Ui::TerrainModelDialog) {
  m_ui->setupUi(this);

  connect(m_ui->generatePushButton, &QPushButton::clicked, this,
          &ccTerrainModelTool::generateAndDisplay);
}

ccTerrainModelTool::~ccTerrainModelTool() { delete m_ui; }

double ccTerrainModelTool::getGridStep() const {
  return m_ui->gridStepDoubleSpinBox->value();
}

unsigned char ccTerrainModelTool::getProjectionDimension() const {
  int dim = m_ui->projDimComboBox->currentIndex();
  assert(dim >= 0 && dim < 3);

  return static_cast<unsigned char>(dim);
}

ccRasterGrid::ProjectionType ccTerrainModelTool::getTypeOfProjection() const {
  switch (m_ui->ProjectionComboBox->currentIndex()) {
  case 0:
    return ccRasterGrid::PROJ_MINIMUM_VALUE;
  case 1:
    return ccRasterGrid::PROJ_MAXIMUM_VALUE;
  case 2:
    return ccRasterGrid::PROJ_MAXIMUM_VALUE;
  default:
    // shouldn't be possible for this option!
    assert(false);
  }

  return ccRasterGrid::INVALID_PROJECTION_TYPE;
}

void ccTerrainModelTool::gridIsUpToDate(bool state) {
  if (state) {
    // standard button
    m_ui->generatePushButton->setStyleSheet(QString());
  } else {
    // red button
    m_ui->generatePushButton->setStyleSheet(
        "color: white; background-color:red;");
  }
  m_ui->generatePushButton->setDisabled(state);
  m_ui->addDBPushButton->setEnabled(state);
  m_ui->exportTiffPushButton->setEnabled(state);
}

void ccTerrainModelTool::saveSettingsAndAccept() {
  saveSettings();
  accept();
}

void ccTerrainModelTool::saveSettings() {
  QSettings settings;
  settings.beginGroup(ccPS::TerrainModel());
  settings.setValue("K", m_ui->spinKNN->value());
  settings.setValue("ProjectionDim", m_ui->projDimComboBox->currentIndex());
  settings.setValue("GridStep", m_ui->gridStepDoubleSpinBox->value());
  settings.setValue("ModelType", m_ui->ProjectionComboBox->currentIndex());
  settings.endGroup();
}

void ccTerrainModelTool::loadSettings() {
  QSettings settings;
  settings.beginGroup(ccPS::VolumeCalculation());
  int knn = settings.value("K", m_ui->spinKNN->value()).toInt();
  int projectionDim =
      settings.value("ProjectionDim", m_ui->projDimComboBox->currentIndex())
          .toInt();
  double step = settings.value("GridStep", m_ui->gridStepDoubleSpinBox->value())
                    .toDouble();
  int modelType =
      settings.value("NumPrecision", m_ui->ProjectionComboBox->currentIndex())
          .toInt();
  settings.endGroup();

  m_ui->spinKNN->setValue(knn);
  m_ui->gridStepDoubleSpinBox->setValue(step);
  m_ui->projDimComboBox->setCurrentIndex(projectionDim);
  m_ui->ProjectionComboBox->setCurrentIndex(modelType);
}

void ccTerrainModelTool::generateAndDisplay() {}

void ccTerrainModelTool::slotGenerate() {
  if (!m_cloud) {
    assert(false);
    return;
  }

  // cloud bounding-box --> grid size
  ccBBox box = getCustomBBox();
  if (!box.isValid()) {
    return;
  }

  unsigned gridWidth = 0;
  unsigned gridHeight = 0;
  if (!getGridSize(gridWidth, gridHeight)) {
    return;
  }

  // grid step
  double gridStep = getGridStep();
  assert(gridStep != 0);

  // cloud
  auto cloud = m_cloud;
  if (cloud == nullptr) {
    assert(false);
    return;
  }

  ccPolyline *edge = nullptr;
  if (m_ui->comboEdge->currentIndex() > 0) {
    edge == m_edge;
  }
}

bool ccTerrainModelTool::generateMode(
    ccRasterGrid &grid, ccGenericPointCloud *ground, ccPolyline *edge,
    const ccBBox &gridBox, unsigned char vertDim, double gridStep,
    unsigned gridWidth, unsigned gridHeight,
    ccRasterGrid::ProjectionType projectionType,
    ccRasterGrid::EmptyCellFillOption groundEmptyCellFillStrategy,
    int krigingKNN, QWidget *parentWidget) {

  return true;
}