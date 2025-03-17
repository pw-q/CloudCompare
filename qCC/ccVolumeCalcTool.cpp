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

#include "ccVolumeCalcTool.h"
#include "ui_volumeCalcDlg.h"

// Local
#include "ccPersistentSettings.h"
#include "mainwindow.h"

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

ccVolumeCalcTool::ccVolumeCalcTool(ccGenericPointCloud *cloud1,
                                   ccGenericPointCloud *cloud2,
                                   ccPolyline *edge,
                                   QWidget *parent /*=nullptr*/)
    : QDialog(parent, Qt::WindowMaximizeButtonHint | Qt::WindowCloseButtonHint),
      cc2Point5DimEditor(), m_cloud1(cloud1), m_cloud2(cloud2), m_edge(edge),
      m_ui(new Ui::VolumeCalcDialog) {
  m_ui->setupUi(this);

  connect(m_ui->buttonBox, &QDialogButtonBox::accepted, this,
          &ccVolumeCalcTool::saveSettingsAndAccept);
  connect(m_ui->buttonBox, &QDialogButtonBox::rejected, this,
          &ccVolumeCalcTool::reject);
  connect(m_ui->gridStepDoubleSpinBox,
          qOverload<double>(&QDoubleSpinBox::valueChanged), this,
          &ccVolumeCalcTool::updateGridInfo);
  connect(m_ui->gridStepDoubleSpinBox,
          qOverload<double>(&QDoubleSpinBox::valueChanged), this,
          &ccVolumeCalcTool::gridOptionChanged);
  connect(m_ui->spinKNN, qOverload<int>(&QSpinBox::valueChanged), this,
          &ccVolumeCalcTool::gridOptionChanged);
  connect(m_ui->groundEmptyValueDoubleSpinBox,
          qOverload<double>(&QDoubleSpinBox::valueChanged), this,
          &ccVolumeCalcTool::gridOptionChanged);
  connect(m_ui->projDimComboBox,
          qOverload<int>(&QComboBox::currentIndexChanged), this,
          &ccVolumeCalcTool::projectionDirChanged);
  connect(m_ui->updatePushButton, &QPushButton::clicked, this,
          &ccVolumeCalcTool::updateGridAndDisplay);
  connect(m_ui->heightProjectionComboBox,
          qOverload<int>(&QComboBox::currentIndexChanged), this,
          &ccVolumeCalcTool::gridOptionChanged);
  connect(m_ui->fillGroundEmptyCellsComboBox,
          qOverload<int>(&QComboBox::currentIndexChanged), this,
          &ccVolumeCalcTool::groundFillEmptyCellStrategyChanged);
  connect(m_ui->swapToolButton, &QToolButton::clicked, this,
          &ccVolumeCalcTool::swapRoles);
  connect(m_ui->groundComboBox, qOverload<int>(&QComboBox::currentIndexChanged),
          this, &ccVolumeCalcTool::groundSourceChanged);
  connect(m_ui->comboEdge, qOverload<int>(&QComboBox::currentIndexChanged),
          this, &ccVolumeCalcTool::edgeChanged);
  connect(m_ui->clipboardPushButton, &QPushButton::clicked, this,
          &ccVolumeCalcTool::exportToClipboard);
  connect(m_ui->exportGridPushButton, &QPushButton::clicked, this,
          &ccVolumeCalcTool::exportGridAsCloud);
  connect(m_ui->precisionSpinBox, qOverload<int>(&QSpinBox::valueChanged), this,
          &ccVolumeCalcTool::setDisplayedNumberPrecision);

  if (m_cloud1 && !m_cloud2) {
    // the existing cloud is always the second by default
    std::swap(m_cloud1, m_cloud2);
  }
  assert(m_cloud2);

  // custom bbox editor
  ccBBox gridBBox = m_cloud1 ? m_cloud1->getOwnBB() : ccBBox();
  if (m_cloud2) {
    gridBBox += m_cloud2->getOwnBB();
  }
  if (gridBBox.isValid()) {
    createBoundingBoxEditor(gridBBox, this);
    connect(m_ui->editGridToolButton, &QToolButton::clicked, this,
            &ccVolumeCalcTool::showGridBoxEditor);
  } else {
    m_ui->editGridToolButton->setEnabled(false);
  }

  m_ui->groundComboBox->addItem("固定值");
  m_ui->ceilComboBox->addItem("固定值");
  m_ui->comboEdge->addItem("None");
  if (m_cloud1) {
    m_ui->groundComboBox->addItem(m_cloud1->getName());
    m_ui->ceilComboBox->addItem(m_cloud1->getName());
  }
  if (m_cloud2) {
    m_ui->groundComboBox->addItem(m_cloud2->getName());
    m_ui->ceilComboBox->addItem(m_cloud2->getName());
  }
  if (m_edge) {
    m_ui->comboEdge->addItem(m_edge->getName());
  }

  m_ui->precisionSpinBox->setVisible(false);
  m_ui->label->setVisible(false);
  m_ui->editGridToolButton->setVisible(false);

  assert(m_ui->groundComboBox->count() >= 2);
  m_ui->groundComboBox->setCurrentIndex(m_ui->groundComboBox->count() - 2);
  m_ui->ceilComboBox->setCurrentIndex(m_ui->ceilComboBox->count() - 1);
  m_ui->comboEdge->setCurrentIndex(m_ui->comboEdge->count() - 1);

  // add window
  create2DView(m_ui->mapFrame);
  if (m_glWindow) {
    ccGui::ParamStruct params = m_glWindow->getDisplayParameters();
    params.colorScaleShowHistogram = false;
    params.displayedNumPrecision = m_ui->precisionSpinBox->value();
    m_glWindow->setDisplayParameters(params, true);
  }

  loadSettings();

  updateGridInfo();

  gridIsUpToDate(false);
}

ccVolumeCalcTool::~ccVolumeCalcTool() { delete m_ui; }

void ccVolumeCalcTool::setDisplayedNumberPrecision(int precision) {
  // update window
  if (m_glWindow) {
    ccGui::ParamStruct params = m_glWindow->getDisplayParameters();
    params.displayedNumPrecision = precision;
    m_glWindow->setDisplayParameters(params, true);
    m_glWindow->redraw(true, false);
  }

  // update report
  if (m_ui->clipboardPushButton->isEnabled()) {
    outputReport(m_lastReport);
  }
}

void ccVolumeCalcTool::groundSourceChanged(int) {
  // m_ui->fillGroundEmptyCellsComboBox->setEnabled(
  //     m_ui->groundComboBox->currentIndex() > 0);
  groundFillEmptyCellStrategyChanged(-1);
}

void ccVolumeCalcTool::edgeChanged(int) { gridIsUpToDate(false); }

void ccVolumeCalcTool::swapRoles() {
  int sourceIndex = m_ui->ceilComboBox->currentIndex();

  m_ui->ceilComboBox->setCurrentIndex(m_ui->groundComboBox->currentIndex());

  m_ui->groundComboBox->setCurrentIndex(sourceIndex);

  gridIsUpToDate(false);
}

bool ccVolumeCalcTool::showGridBoxEditor() {
  if (cc2Point5DimEditor::showGridBoxEditor()) {
    updateGridInfo();
    return true;
  }

  return false;
}

void ccVolumeCalcTool::gridOptionChanged() { gridIsUpToDate(false); }

void ccVolumeCalcTool::updateGridInfo() {
  m_ui->gridWidthLabel->setText(getGridSizeAsString());
}

double ccVolumeCalcTool::getGridStep() const {
  return m_ui->gridStepDoubleSpinBox->value();
}

unsigned char ccVolumeCalcTool::getProjectionDimension() const {
  int dim = m_ui->projDimComboBox->currentIndex();
  assert(dim >= 0 && dim < 3);

  return static_cast<unsigned char>(dim);
}

void ccVolumeCalcTool::sfProjectionTypeChanged(int index) {
  Q_UNUSED(index)

  gridIsUpToDate(false);
}

void ccVolumeCalcTool::projectionDirChanged(int dir) {
  Q_UNUSED(dir)

  updateGridInfo();
  gridIsUpToDate(false);
}

void ccVolumeCalcTool::groundFillEmptyCellStrategyChanged(int) {
  ccRasterGrid::EmptyCellFillOption fillEmptyCellsStrategy =
      getFillEmptyCellsStrategy(m_ui->fillGroundEmptyCellsComboBox);

  m_ui->groundEmptyValueDoubleSpinBox->setEnabled(
      (m_ui->groundComboBox->currentIndex() == 0) ||
      (fillEmptyCellsStrategy == ccRasterGrid::FILL_CUSTOM_HEIGHT));

  m_ui->spinKNN->setEnabled(fillEmptyCellsStrategy ==
                            ccRasterGrid::INTERPOLATE_KRIGING);

  gridIsUpToDate(false);
}

ccRasterGrid::ProjectionType ccVolumeCalcTool::getTypeOfProjection() const {
  switch (m_ui->heightProjectionComboBox->currentIndex()) {
  case 0:
    return ccRasterGrid::PROJ_MINIMUM_VALUE;
  case 1:
    return ccRasterGrid::PROJ_AVERAGE_VALUE;
  case 2:
    return ccRasterGrid::PROJ_MAXIMUM_VALUE;
  default:
    // shouldn't be possible for this option!
    assert(false);
  }

  return ccRasterGrid::INVALID_PROJECTION_TYPE;
}

void ccVolumeCalcTool::loadSettings() {
  QSettings settings;
  settings.beginGroup(ccPS::VolumeCalculation());
  int projType = settings
                     .value("ProjectionType",
                            m_ui->heightProjectionComboBox->currentIndex())
                     .toInt();
  int projDim =
      settings.value("ProjectionDim", m_ui->projDimComboBox->currentIndex())
          .toInt();
  int groundFillStrategy =
      settings
          .value("gFillStrategy",
                 m_ui->fillGroundEmptyCellsComboBox->currentIndex())
          .toInt();
  double step = settings.value("GridStep", m_ui->gridStepDoubleSpinBox->value())
                    .toDouble();
  double groundEmptyHeight =
      settings
          .value("gEmptyCellsHeight",
                 m_ui->groundEmptyValueDoubleSpinBox->value())
          .toDouble();
  double groundMaxEdgeLength =
      settings.value("gMaxEdgeLength", m_ui->spinKNN->value()).toDouble();
  int precision =
      settings.value("NumPrecision", m_ui->precisionSpinBox->value()).toInt();
  settings.endGroup();

  m_ui->gridStepDoubleSpinBox->setValue(step);
  m_ui->heightProjectionComboBox->setCurrentIndex(projType);
  m_ui->fillGroundEmptyCellsComboBox->setCurrentIndex(groundFillStrategy);
  m_ui->groundEmptyValueDoubleSpinBox->setValue(groundEmptyHeight);
  m_ui->spinKNN->setValue(groundMaxEdgeLength);
  m_ui->projDimComboBox->setCurrentIndex(projDim);
  m_ui->precisionSpinBox->setValue(precision);
}

void ccVolumeCalcTool::saveSettingsAndAccept() {
  saveSettings();
  accept();
}

void ccVolumeCalcTool::saveSettings() {
  QSettings settings;
  settings.beginGroup(ccPS::VolumeCalculation());
  settings.setValue("ProjectionType",
                    m_ui->heightProjectionComboBox->currentIndex());
  settings.setValue("ProjectionDim", m_ui->projDimComboBox->currentIndex());
  settings.setValue("gFillStrategy",
                    m_ui->fillGroundEmptyCellsComboBox->currentIndex());
  settings.setValue("GridStep", m_ui->gridStepDoubleSpinBox->value());
  settings.setValue("gEmptyCellsHeight",
                    m_ui->groundEmptyValueDoubleSpinBox->value());
  settings.setValue("gMaxEdgeLength", m_ui->spinKNN->value());
  settings.setValue("NumPrecision", m_ui->precisionSpinBox->value());
  settings.endGroup();
}

void ccVolumeCalcTool::gridIsUpToDate(bool state) {
  if (state) {
    // standard button
    m_ui->updatePushButton->setStyleSheet(QString());
  } else {
    // red button
    m_ui->updatePushButton->setStyleSheet(
        "color: white; background-color:red;");
  }
  m_ui->updatePushButton->setDisabled(state);
  m_ui->clipboardPushButton->setEnabled(state);
  m_ui->exportGridPushButton->setEnabled(state);
  if (!state) {
    m_ui->spareseWarningLabel->hide();
    m_ui->reportPlainTextEdit->setPlainText("请先设置参数，进行体积计算");
  }
}

ccPointCloud *ccVolumeCalcTool::ConvertGridToCloud(ccRasterGrid &grid,
                                                   const ccBBox &gridBox,
                                                   unsigned char vertDim,
                                                   bool exportToOriginalCS) {
  assert(gridBox.isValid());
  assert(vertDim < 3);

  ccPointCloud *rasterCloud = nullptr;
  try {
    // we only compute the default 'height' layer
    std::vector<ccRasterGrid::ExportableFields> exportedStatistics(1);
    exportedStatistics.front() = ccRasterGrid::PER_CELL_VALUE;

    rasterCloud = grid.convertToCloud(true, false, exportedStatistics, false,
                                      false, false, false, nullptr, vertDim,
                                      gridBox, 0, exportToOriginalCS, false);

    if (rasterCloud && rasterCloud->hasScalarFields()) {
      rasterCloud->showSF(true);
      rasterCloud->setCurrentDisplayedScalarField(0);
      ccScalarField *sf =
          static_cast<ccScalarField *>(rasterCloud->getScalarField(0));
      assert(sf);
      sf->setName("Relative height");
      sf->setSymmetricalScale(sf->getMin() < 0 && sf->getMax() > 0);
      rasterCloud->showSFColorsScale(true);
    }
  } catch (const std::bad_alloc &) {
    ccLog::Warning("[ConvertGridToCloud] 没有足够的内存!");
    if (rasterCloud) {
      delete rasterCloud;
      rasterCloud = nullptr;
    }
  }

  return rasterCloud;
}

ccPointCloud *
ccVolumeCalcTool::convertGridToCloud(bool exportToOriginalCS) const {
  ccPointCloud *rasterCloud = nullptr;
  try {
    // we only compute the default 'height' layer
    std::vector<ccRasterGrid::ExportableFields> exportedStatistics(1);
    exportedStatistics.front() = ccRasterGrid::PER_CELL_VALUE;
    rasterCloud = cc2Point5DimEditor::convertGridToCloud(
        true, false, exportedStatistics, false, false, false, false, nullptr,
        0.0, exportToOriginalCS, false, nullptr);

    if (rasterCloud) {
      if (rasterCloud->hasScalarFields()) {
        rasterCloud->showSF(true);
        rasterCloud->setCurrentDisplayedScalarField(0);
        ccScalarField *sf =
            static_cast<ccScalarField *>(rasterCloud->getScalarField(0));
        assert(sf);
        sf->setName("Relative height");
        sf->setSymmetricalScale(sf->getMin() < 0 && sf->getMax() > 0);
        rasterCloud->showSFColorsScale(true);
      }

      // keep Global Shift & Scale
      auto ground = getGroundCloud();
      auto ceil = getCeilCloud();

      if (ground.first && ground.first->isShifted()) {
        rasterCloud->copyGlobalShiftAndScale(*ground.first);
      } else if (ceil.first && ceil.first->isShifted()) {
        rasterCloud->copyGlobalShiftAndScale(*ceil.first);
      }
    }
  } catch (const std::bad_alloc &) {
    ccLog::Error("没有足够内存!");
    if (rasterCloud) {
      delete rasterCloud;
      rasterCloud = nullptr;
    }
  }

  return rasterCloud;
}

void ccVolumeCalcTool::updateGridAndDisplay() {
  bool success = updateGrid();
  if (success && m_glWindow) {
    // convert grid to point cloud
    if (m_rasterCloud) {
      m_glWindow->removeFromOwnDB(m_rasterCloud);
      delete m_rasterCloud;
      m_rasterCloud = nullptr;
    }
    if (m_rasterLine) {
      m_glWindow->removeFromOwnDB(m_rasterLine);
      delete m_rasterLine;
      m_rasterLine = nullptr;
    }

    m_rasterCloud = convertGridToCloud(false);
    if (m_rasterCloud) {
      m_glWindow->addToOwnDB(m_rasterCloud);
      ccBBox box = m_rasterCloud->getDisplayBB_recursive(false, m_glWindow);
      update2DDisplayZoom(box);
    } else {
      ccLog::Error("没有足够内存!");
      m_glWindow->redraw();
    }
    if (m_edge) {
      m_rasterLine = new ccPolyline(*m_edge);
      m_rasterLine->setWidth(2);
      m_rasterLine->setColor(ccColor::red);
      m_rasterLine->showColors(true);
      m_rasterLine->set2DMode(false);
      m_rasterLine->setSelected(false);

      m_glWindow->addToOwnDB(m_rasterLine);
    } //
  }

  gridIsUpToDate(success);
}

QString ccVolumeCalcTool::ReportInfo::toText(int precision) const {
  QLocale locale(QLocale::English);

  QStringList reportText;
  reportText
      << QString("体积: %1").arg(locale.toString(volume, 'f', precision));
  reportText
      << QString("表面积: %1").arg(locale.toString(surface, 'f', precision));
  reportText << QString("----------------------");
  reportText << QString("增加的体积: (+)%1")
                    .arg(locale.toString(addedVolume, 'f', precision));
  reportText << QString("移除的体积: (-)%1")
                    .arg(locale.toString(removedVolume, 'f', precision));
  reportText << QString("----------------------");
  reportText << QString("匹配格网: %1%").arg(matchingPrecent, 0, 'f', 1);
  reportText << QString("不匹配格网:");
  reportText
      << QString("    before: %1%").arg(groundNonMatchingPercent, 0, 'f', 1);
  reportText
      << QString("    after: %1%").arg(ceilNonMatchingPercent, 0, 'f', 1);
  reportText << QString("平均邻居数量: %1 / 8.0")
                    .arg(averageNeighborsPerCell, 0, 'f', 1);

  return reportText.join("\n");
}

void ccVolumeCalcTool::outputReport(const ReportInfo &info) {
  int precision = m_ui->precisionSpinBox->value();

  m_ui->reportPlainTextEdit->setPlainText(info.toText(precision));

  // below 7 neighbors per cell, at least one of the cloud is very sparse!
  m_ui->spareseWarningLabel->setVisible(info.averageNeighborsPerCell < 7.0f);

  m_lastReport = info;
  m_ui->clipboardPushButton->setEnabled(true);
}

bool SendError(const QString &message, QWidget *parentWidget) {
  if (parentWidget) {
    ccLog::Error(message);
  } else {
    ccLog::Warning("[体积计算] " + message);
  }
  return false;
}

bool ccVolumeCalcTool::ComputeVolume(
    ccRasterGrid &grid, ccGenericPointCloud *ground, ccGenericPointCloud *ceil,
    ccPolyline *edge, const ccBBox &gridBox, unsigned char vertDim,
    double gridStep, unsigned gridWidth, unsigned gridHeight,
    ccRasterGrid::ProjectionType projectionType,
    ccRasterGrid::EmptyCellFillOption groundEmptyCellFillStrategy,
    int krigingKNN, ccRasterGrid::EmptyCellFillOption ceilEmptyCellFillStrategy,
    double ceilMaxEdgeLength, ccVolumeCalcTool::ReportInfo &reportInfo,
    double groundHeight = std::numeric_limits<double>::quiet_NaN(),
    double ceilHeight = std::numeric_limits<double>::quiet_NaN(),
    QWidget *parentWidget /*=nullptr*/) {
  if (gridStep <= 1.0e-8 || gridWidth == 0 || gridHeight == 0 || vertDim > 2) {
    assert(false);
    ccLog::Warning("[体积计算] 无效输入参数");
    return false;
  }

  if (!ground && !ceil) {
    assert(false);
    ccLog::Warning("[体积计算] 无效输入点云");
    return false;
  }

  if (!gridBox.isValid()) {
    ccLog::Warning("[体积计算] 无效包围盒");
    return false;
  }

  // grid size
  unsigned gridTotalSize = gridWidth * gridHeight;
  if (gridTotalSize == 1) {
    if (parentWidget &&
        QMessageBox::question(parentWidget, "格网大小异常",
                              "格网大小仅为1"
                              "是否继续处理?",
                              QMessageBox::Yes,
                              QMessageBox::No) == QMessageBox::No)
      return false;
  } else if (gridTotalSize > 10000000) {
    if (parentWidget &&
        QMessageBox::question(parentWidget, "格网大小异常",
                              "生成格网数量已大于10000000, "
                              "是否继续处理?",
                              QMessageBox::Yes,
                              QMessageBox::No) == QMessageBox::No)
      return false;
  }

  // memory allocation
  CCVector3d minCorner = gridBox.minCorner();
  if (!grid.init(gridWidth, gridHeight, gridStep, minCorner)) {
    // not enough memory
    return SendError("没有足够内存", parentWidget);
  }

  // progress dialog
  QScopedPointer<ccProgressDialog> pDlg(nullptr);
  if (parentWidget) {
    pDlg.reset(new ccProgressDialog(true, parentWidget));
  }

  ccRasterGrid groundRaster;
  if (ground) {
    if (!groundRaster.init(gridWidth, gridHeight, gridStep, minCorner)) {
      // not enough memory
      return SendError("没有足够内存", parentWidget);
    }

    ccRasterGrid::InterpolationType interpolationType =
        ccRasterGrid::InterpolationTypeFromEmptyCellFillOption(
            groundEmptyCellFillStrategy);
    // ccRasterGrid::DelaunayInterpolationParams dInterpParams;
    ccRasterGrid::KrigingParams dKrigingParams;
    void *interpolationParams = nullptr;
    switch (interpolationType) {
    case ccRasterGrid::InterpolationType::DELAUNAY:
      // dInterpParams.maxEdgeLength = groundMaxEdgeLength;
      // interpolationParams = (void *)&dInterpParams;
      break;
    case ccRasterGrid::InterpolationType::KRIGING:
      dKrigingParams.kNN = krigingKNN;
      interpolationParams = (void *)&dKrigingParams;
      // not supported yet
      // assert(false);
      break;
    default:
      // do nothing
      break;
    }

    if (groundRaster.fillWith(ground, vertDim, projectionType,
                              interpolationType, interpolationParams,
                              ccRasterGrid::INVALID_PROJECTION_TYPE,
                              pDlg.data())) {
      groundRaster.fillEmptyCells(groundEmptyCellFillStrategy, groundHeight);
      ccLog::Print(
          QString(
              "[体积计算] 参考点云格网信息: 大小: %1 x %2 / 高度: [%3 ; %4]")
              .arg(groundRaster.width)
              .arg(groundRaster.height)
              .arg(groundRaster.minHeight)
              .arg(groundRaster.maxHeight));
    } else {
      return false;
    }
  }

  // ceil
  ccRasterGrid ceilRaster;
  if (ceil) {
    if (!ceilRaster.init(gridWidth, gridHeight, gridStep, minCorner)) {
      // not enough memory
      return SendError("没有足够内存", parentWidget);
    }

    ccRasterGrid::InterpolationType interpolationType =
        ccRasterGrid::InterpolationTypeFromEmptyCellFillOption(
            ceilEmptyCellFillStrategy);
    ccRasterGrid::DelaunayInterpolationParams dInterpParams;
    ccRasterGrid::KrigingParams dKrigingParams;
    void *interpolationParams = nullptr;
    switch (interpolationType) {
    case ccRasterGrid::InterpolationType::DELAUNAY:
      dInterpParams.maxEdgeLength = ceilMaxEdgeLength;
      interpolationParams = (void *)&dInterpParams;
      break;
    case ccRasterGrid::InterpolationType::KRIGING:
      dKrigingParams.kNN = krigingKNN;
      interpolationParams = (void *)&dKrigingParams;
      // not supported yet
      // assert(false);
      break;
    default:
      // do nothing
      break;
    }

    if (ceilRaster.fillWith(ceil, vertDim, projectionType, interpolationType,
                            interpolationParams,
                            ccRasterGrid::INVALID_PROJECTION_TYPE,
                            pDlg.data())) {
      ceilRaster.fillEmptyCells(ceilEmptyCellFillStrategy, ceilHeight);
      ccLog::Print(
          QString(
              "[体积计算] 计算点云格网信息: 大小: %1 x %2 / 高度: [%3 ; %4]")
              .arg(ceilRaster.width)
              .arg(ceilRaster.height)
              .arg(ceilRaster.minHeight)
              .arg(ceilRaster.maxHeight));
    } else {
      return false;
    }
  }

  // update grid and compute volume
  {
    if (pDlg) {
      pDlg->setMethodTitle(QObject::tr("体积计算"));
      pDlg->setInfo(
          QObject::tr("格网: %1 x %2").arg(grid.width).arg(grid.height));
      pDlg->start();
      pDlg->show();
      QCoreApplication::processEvents();
    }
    CCCoreLib::NormalizedProgress nProgress(pDlg.data(),
                                            grid.width * grid.height);

    size_t ceilNonMatchingCount = 0;
    size_t groundNonMatchingCount = 0;
    size_t cellCount = 0;

    if (edge != nullptr)
      grid.getInterRatioWithPolyline(edge);

    // at least one of the grid is based on a cloud
    grid.nonEmptyCellCount = 0;
    for (unsigned i = 0; i < grid.height; ++i) {
      for (unsigned j = 0; j < grid.width; ++j) {
        ccRasterCell &cell = grid.rows[i][j];

        bool validGround = true;
        cell.minHeight = groundHeight;
        if (ground) {
          cell.minHeight = groundRaster.rows[i][j].h;
          validGround =
              std::isfinite(cell.minHeight) && cell.interAreaRatio != 0;
        }

        bool validCeil = true;
        cell.maxHeight = ceilHeight;
        if (ceil) {
          cell.maxHeight = ceilRaster.rows[i][j].h;
          validCeil =
              (std::isfinite(cell.maxHeight) && cell.interAreaRatio != 0);
        }

        if (validGround && validCeil) {
          cell.h = cell.maxHeight - cell.minHeight;
          cell.nbPoints = 1;

          cell.h *= cell.interAreaRatio;
          reportInfo.volume += cell.h;
          if (cell.h < 0) {
            reportInfo.removedVolume -= cell.h;
          } else if (cell.h > 0) {
            reportInfo.addedVolume += cell.h;
          }
          reportInfo.surface += 1.0 * cell.interAreaRatio;
          ++grid.nonEmptyCellCount; // matching count
          ++cellCount;
        } else {
          if (validGround) {
            ++cellCount;
            ++groundNonMatchingCount;
          } else if (validCeil) {
            ++cellCount;
            ++ceilNonMatchingCount;
          }
          cell.h = std::numeric_limits<double>::quiet_NaN();
          cell.nbPoints = 0;
        }

        if (pDlg && !nProgress.oneStep()) {
          ccLog::Warning("[体积计算]用户取消计算");
          return false;
        }
      }
    }
    grid.validCellCount = grid.nonEmptyCellCount;

    // count the average number of valid neighbors
    {
      size_t validNeighborsCount = 0;
      size_t count = 0;
      for (unsigned i = 1; i < grid.height - 1; ++i) {
        for (unsigned j = 1; j < grid.width - 1; ++j) {
          ccRasterCell &cell = grid.rows[i][j];
          if (std::isfinite(cell.h)) {
            for (unsigned k = i - 1; k <= i + 1; ++k) {
              for (unsigned l = j - 1; l <= j + 1; ++l) {
                if (k != i || l != j) {
                  ccRasterCell &otherCell = grid.rows[k][l];
                  if (std::isfinite(otherCell.h)) {
                    ++validNeighborsCount;
                  }
                }
              }
            }

            ++count;
          }
        }
      }

      if (count) {
        reportInfo.averageNeighborsPerCell =
            static_cast<double>(validNeighborsCount) / count;
      }
    }

    reportInfo.matchingPrecent =
        static_cast<float>(grid.validCellCount * 100) / cellCount;
    reportInfo.groundNonMatchingPercent =
        static_cast<float>(groundNonMatchingCount * 100) / cellCount;
    reportInfo.ceilNonMatchingPercent =
        static_cast<float>(ceilNonMatchingCount * 100) / cellCount;
    float cellArea = static_cast<float>(grid.gridStep * grid.gridStep);
    reportInfo.volume *= cellArea;
    reportInfo.addedVolume *= cellArea;
    reportInfo.removedVolume *= cellArea;
    reportInfo.surface *= cellArea;
  }

  grid.setValid(true);

  return true;
}

std::pair<ccGenericPointCloud *, double>
ccVolumeCalcTool::getGroundCloud() const {
  ccGenericPointCloud *groundCloud = nullptr;
  double groundHeight = std::numeric_limits<double>::quiet_NaN();
  switch (m_ui->groundComboBox->currentIndex()) {
  case 0:
    groundHeight = m_ui->groundEmptyValueDoubleSpinBox->value();
    break;
  case 1:
    groundCloud = m_cloud1 ? m_cloud1 : m_cloud2;
    break;
  case 2:
    groundCloud = m_cloud2;
    break;
  default:
    assert(false);
    break;
  }

  return {groundCloud, groundHeight};
}

std::pair<ccGenericPointCloud *, double>
ccVolumeCalcTool::getCeilCloud() const {
  ccGenericPointCloud *ceilCloud = nullptr;
  double ceilHeight = std::numeric_limits<double>::quiet_NaN();
  switch (m_ui->ceilComboBox->currentIndex()) {
  case 0:
    ceilHeight = m_ui->groundEmptyValueDoubleSpinBox->value();
    break;
  case 1:
    ceilCloud = m_cloud1 ? m_cloud1 : m_cloud2;
    break;
  case 2:
    ceilCloud = m_cloud2;
    break;
  default:
    assert(false);
    break;
  }

  return {ceilCloud, ceilHeight};
}

bool ccVolumeCalcTool::updateGrid() {
  if (!m_cloud2) {
    assert(false);
    return false;
  }

  // cloud bounding-box --> grid size
  ccBBox box = getCustomBBox();
  if (!box.isValid()) {
    return false;
  }

  unsigned gridWidth = 0;
  unsigned gridHeight = 0;
  if (!getGridSize(gridWidth, gridHeight)) {
    return false;
  }

  // grid step
  double gridStep = getGridStep();
  assert(gridStep != 0);

  // ground
  auto ground = getGroundCloud();
  if (nullptr == ground.first && std::isnan(ground.second)) {
    assert(false);
    return false;
  }

  // ceil
  auto ceil = getCeilCloud();
  if (nullptr == ceil.first && std::isnan(ceil.second)) {
    assert(false);
    return false;
  }

  ccPolyline *edge = nullptr;
  if (m_ui->comboEdge->currentIndex() > 0) {
    edge == m_edge;
  }

  ccVolumeCalcTool::ReportInfo reportInfo;

  if (ComputeVolume(
          m_grid, ground.first, ceil.first, edge, box, getProjectionDimension(),
          gridStep, gridWidth, gridHeight, getTypeOfProjection(),
          getFillEmptyCellsStrategy(m_ui->fillGroundEmptyCellsComboBox),
          m_ui->spinKNN->value(),
          getFillEmptyCellsStrategy(m_ui->fillGroundEmptyCellsComboBox),
          m_ui->spinKNN->value(), reportInfo, ground.second, ceil.second,
          this)) {
    outputReport(reportInfo);
    return true;
  } else {
    return false;
  }
}

void ccVolumeCalcTool::exportToClipboard() const {
  QClipboard *clipboard = QApplication::clipboard();
  if (clipboard) {
    clipboard->setText(m_ui->reportPlainTextEdit->toPlainText());
  }
}

void ccVolumeCalcTool::exportGridAsCloud() const {
  if (!m_grid.isValid()) {
    assert(false);
  }

  ccPointCloud *rasterCloud = convertGridToCloud(true);
  if (!rasterCloud) {
    // error message should have already been issued
    return;
  }

  rasterCloud->setName("高度差 " + rasterCloud->getName());
  ccGenericPointCloud *originCloud = (m_cloud1 ? m_cloud1 : m_cloud2);
  assert(originCloud);
  if (originCloud) {
    if (originCloud->getParent()) {
      originCloud->getParent()->addChild(rasterCloud);
    }
    rasterCloud->setDisplay(originCloud->getDisplay());
  }

  MainWindow *mainWindow = MainWindow::TheInstance();
  if (mainWindow) {
    mainWindow->addToDB(rasterCloud);
    ccLog::Print(
        QString("[Volume] 点云 '%1' 成功导出").arg(rasterCloud->getName()));
  } else {
    assert(false);
    delete rasterCloud;
  }
}
