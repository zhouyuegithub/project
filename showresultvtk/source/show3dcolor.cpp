/**********************************************************************

  文件名: 6.7_PolyDataMarchingCubes.cpp
  Copyright (c) 张晓东, 罗火灵. All rights reserved.
  更多信息请访问: 
    http://www.vtkchina.org (VTK中国)
	http://blog.csdn.net/www_doling_net (东灵工作室) 

**********************************************************************/

#include <vtkSmartPointer.h>
#include <vtkMarchingCubes.h>
#include <vtkContourFilter.h>
#include <vtkVoxelModeller.h>
#include <vtkSphereSource.h>
#include <vtkImageData.h>
#include <vtkMetaImageReader.h>
#include <vtkInteractorStyleImage.h>
#include <vtkActor.h>
#include <vtkPolyDataMapper.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkImageActor.h>
#include <vtkProperty.h>
#include <cassert>
#include <string>

#include <vector>
#include "itksys/SystemTools.hxx"

#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "itkBSplineScatteredDataPointSetToImageFilter.h"
#include "p.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageToVTKImageFilter.h"
#include "itkGDCMImageIO.h"
#include "itkGDCMSeriesFileNames.h"
#include "itkImageSeriesReader.h"
#include <fstream>
using namespace std;

#include "itkGDCMImageIO.h"
#include "itkGDCMSeriesFileNames.h"
#include "itkImageSeriesReader.h"
#include "itkImageSeriesWriter.h"
#include "itkCommand.h"
#include "itkVersorRigid3DTransform.h"
#include "itkCenteredTransformInitializer.h"
#include "itkRegularStepGradientDescentOptimizerv4.h"
#include "itkThresholdImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkHistogramMatchingImageFilter.h"
#include "itkImageFileWriter.h"
#include "itkFEMRegistrationFilter.h"
#include "itkMeanSquaresImageToImageMetricv4.h"
#include "itkImageRegistrationMethodv4.h"
#include "itkRegularStepGradientDescentOptimizerv4.h"
#include "itkVersorRigid3DTransformOptimizer.h"
#include "itkVersorRigid3DTransform.h"
#include "itkCenteredTransformInitializer.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkResampleImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkSubtractImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkExtractImageFilter.h"
//测试：../data/HeadMRVolume.mhd 200
const char* confixpth = "C:\\Users\\Administrator\\Desktop\\devide\\zhangmengmeng500\\concon\\8zhangmengmeng500SEGandFILL";
const char* conmovepth = "C:\\Users\\Administrator\\Desktop\\devide\\zhangmengmeng500\\concon\\2zhangnmengmeng500SEGandFILL";
typedef itk::VersorRigid3DTransform< double > TransformType;
TransformType::Pointer finalTransform = TransformType::New();
TransformType::MatrixType matrix ;
TransformType::OffsetType offset ;
int main(/*int argc, char *argv[]*/)
{
	ifstream fin;
	string txtptn = "C:\\Users\\Administrator\\Desktop\\result\\zhangmengmeng\\1st\\2.txt";
	fin.open(txtptn.data());
	//TransformType::MatrixType matrix ;
	//TransformType::OffsetType offset ;
	string s;
	while (getline(fin,s))
	{
		if(s=="Offset1")
		{
			getline(fin,s);//nextline
			//djksdjkh
			finalTransform->SetOffset(offset);

		}
	}
	typedef signed short    Pixel;	
	typedef itk::Image< Pixel, 3 >         ImageType;
	typedef itk::ImageSeriesReader< ImageType >        ReaderType;
	ReaderType::Pointer T1reader = ReaderType::New();//定义一个T1reader的指针用于读图
	typedef itk::GDCMImageIO       ImageIOType;
	ImageIOType::Pointer T1dicomIO = ImageIOType::New();///定义T1dicomIO指针
	T1reader->SetImageIO( T1dicomIO );//实际执行读取任务的ImageIO对象现在被连接到ImageSeriesReader。这是确保我们用适于我们想要读取的文件类型的ImageIO对象的最安全方法
	typedef itk::GDCMSeriesFileNames NamesGeneratorType;
	NamesGeneratorType::Pointer T1nameGenerator = NamesGeneratorType::New();
	T1nameGenerator->SetUseSeriesDetails( true );
	T1nameGenerator->SetDirectory(conmovepth);
	typedef std::vector< std::string >   FileNamesContainer;
	try
	{
		std::cout << std::endl << "The directory: " << std::endl;
		std::cout << std::endl << "first"<< std::endl << std::endl;
		std::cout << "Contains the following DICOM Series: ";
		std::cout << std::endl << std::endl;
		typedef std::vector< std::string >    SeriesIdContainer;
		//SeriesIdContainer类型的T1seriesUID存储序列信息用于唯一标识DICOM标准中各种不同信息对象
		const SeriesIdContainer & T1seriesUID = T1nameGenerator->GetSeriesUIDs();
		SeriesIdContainer::const_iterator T1seriesItr = T1seriesUID.begin();////迭代器
		SeriesIdContainer::const_iterator T1seriesEnd = T1seriesUID.end();
		//显示读取到的UID
		while( T1seriesItr != T1seriesEnd )
		{
			std::cout << T1seriesItr->c_str() << std::endl;
			++T1seriesItr;
		}
		std::string T1seriesIdentifier;

		T1seriesIdentifier = T1seriesUID.begin()->c_str();//通过迭代器读取所有单张切片
		std::cout << std::endl << std::endl;
		std::cout << "Now reading series: " << std::endl << std::endl;
		std::cout << T1seriesIdentifier << std::endl;
		std::cout << std::endl << std::endl;
		FileNamesContainer T1fileNames;
		T1fileNames = T1nameGenerator->GetFileNames( T1seriesIdentifier );///读取文件
		T1reader->SetFileNames( T1fileNames );
		//最后，我们用writer的Update( )触发管道的执行。这时图像的切片保存在单独的文件中，每个文件包含一个单独的切片。用于这些切片的文件名由文件名发生器生成
		try
		{
			T1reader->Update();
		}
		catch (itk::ExceptionObject &ex)
		{
			std::cout << ex << std::endl;
			return EXIT_FAILURE;
		}
	}
	catch (itk::ExceptionObject &ex)
	{
		std::cout << ex << std::endl;
		return EXIT_FAILURE;
	}

	////////////////////////////////使用itk加载第二幅序列//////////////////////
	ReaderType::Pointer T2reader = ReaderType::New();
	ImageIOType::Pointer T2dicomIO = ImageIOType::New();
	T2reader->SetImageIO( T2dicomIO );
	NamesGeneratorType::Pointer T2nameGenerator = NamesGeneratorType::New();//读取时的文件名有文件名发生器产生
	T2nameGenerator->SetUseSeriesDetails( true );
//	T2nameGenerator->AddSeriesRestriction("0008|0021" );
	/*	T2nameGenerator->SetDirectory( "D:\\anonymized_before\\0020\\3" );*/
	//T2nameGenerator->SetDirectory( "C:\\Users\\Administrator\\Desktop\\devide\\haolinlin500\\condyle\\12haolinlinSEGandFILL"/*"C:\\Users\\Administrator\\Desktop\\devide\\square2" */);
	T2nameGenerator->SetDirectory(confixpth);
	//T2nameGenerator->SetDirectory( "C:\\Users\\Administrator\\Desktop\\niuxiahe\\mha文件读取\\新建文件夹\\2" );
	try
	{
		std::cout << std::endl << "The directory: " << std::endl;
		std::cout << std::endl << "16chengSEG"<< std::endl << std::endl;
		std::cout << "Contains the following DICOM Series: ";///每个dicom序列的信息有UID提供
		std::cout << std::endl << std::endl;
		typedef std::vector< std::string >    SeriesIdContainer;
		const SeriesIdContainer & T2seriesUID = T2nameGenerator->GetSeriesUIDs();////????????????????????????
		SeriesIdContainer::const_iterator T2seriesItr = T2seriesUID.begin();
		SeriesIdContainer::const_iterator T2seriesEnd = T2seriesUID.end();
		while( T2seriesItr != T2seriesEnd )
		{
			std::cout << T2seriesItr->c_str() << std::endl;
			++T2seriesItr;
		}
		std::string T2seriesIdentifier;

		T2seriesIdentifier = T2seriesUID.begin()->c_str();

		std::cout << std::endl << std::endl;
		std::cout << "Now reading series: " << std::endl << std::endl;
		std::cout << T2seriesIdentifier << std::endl;
		std::cout << std::endl << std::endl;
		FileNamesContainer T2fileNames;
		T2fileNames = T2nameGenerator->GetFileNames( T2seriesIdentifier );
		T2reader->SetFileNames( T2fileNames );
		try
		{
			T2reader->Update();
		}
		catch (itk::ExceptionObject &ex)
		{
			std::cout << ex << std::endl;
			return EXIT_FAILURE;
		}
	}
	catch (itk::ExceptionObject &ex)
	{
		std::cout << ex << std::endl;
		return EXIT_FAILURE;
	}
	const unsigned int      Dimension = 3;
	typedef  float            FloatPixelType;
	typedef itk::Image< FloatPixelType, Dimension >      FloatImageType;////输入已经定义过，这是输出
	typedef itk::CastImageFilter<
	ImageType, FloatImageType >  CastFilterType;
//CastFilterType::Pointer       allT2_castFilter       = CastFilterType::New();
CastFilterType::Pointer       DWI_castFilter       = CastFilterType::New();
CastFilterType::Pointer       T1_castFilter       = CastFilterType::New();
		T1_castFilter->SetInput(       T1reader->GetOutput() );///将读取到的数据T1reader通过castfiltertype得到滤波后的图像
	T1_castFilter->Update();

	//CastFilterType::Pointer       DWI_castFilter       = CastFilterType::New();
	DWI_castFilter->SetInput(       T2reader->GetOutput() );///将read的图作为输入
	DWI_castFilter->Update();

	/*typedef itk::RescaleIntensityImageFilter<
	FloatImageType, FloatImageType > FloatRescaleFilterType;



	FloatRescaleFilterType::Pointer T1_rescalercast = FloatRescaleFilterType::New();*/
	typedef itk::RescaleIntensityImageFilter<
	FloatImageType, FloatImageType > FloatRescaleFilterType;


FloatRescaleFilterType::Pointer DWI_filter = FloatRescaleFilterType::New();
	FloatRescaleFilterType::Pointer T1_filter = FloatRescaleFilterType::New();
	T1_filter->SetOutputMinimum(   0 );
	T1_filter->SetOutputMaximum( 255 );
	T1_filter->SetInput( T1_castFilter->GetOutput());
	T1_filter->Update();
	//FloatRescaleFilterType::Pointer DWI_rescalercast = FloatRescaleFilterType::New();
	DWI_filter->SetOutputMinimum(   0 );
	DWI_filter->SetOutputMaximum( 255 );///归到0-255
	DWI_filter->SetInput( DWI_castFilter->GetOutput());
	DWI_filter->Update();
	typedef itk::ImageToVTKImageFilter< FloatImageType> itkTovtkFilterType;
	itkTovtkFilterType::Pointer itkTovtkImageFilter1 = itkTovtkFilterType::New();
	itkTovtkImageFilter1->SetInput(T1_castFilter->GetOutput());
	itkTovtkImageFilter1->Update();

	double isoValue = 255;

	vtkSmartPointer<vtkMarchingCubes> surface =
	    vtkSmartPointer<vtkMarchingCubes>::New();//对数据进行MarchingCubes 算法的处理，输出重建的物体表面数据。
	surface->SetInputData(itkTovtkImageFilter1->GetOutput());
	surface->ComputeNormalsOn();//计算等值面法向量提高渲染质量
	surface->SetValue(0, isoValue); 
  
	vtkSmartPointer<vtkPolyDataMapper> surfMapper =
	    vtkSmartPointer<vtkPolyDataMapper>::New();
	surfMapper->SetInputConnection(surface->GetOutputPort());//渲染多边形几何数据
	surfMapper->ScalarVisibilityOff();

	vtkSmartPointer<vtkActor> surfActor =
	    vtkSmartPointer<vtkActor>::New();//渲染场景中数据的可视化表达显示图像
	surfActor->SetMapper(surfMapper);
	surfActor->GetProperty()->SetColor(0.0, 0.0, 1.0);

	vtkSmartPointer<vtkRenderer> renderer =
	   vtkSmartPointer<vtkRenderer>::New();//负责管理场景的渲染过程将场景对象组合起来
	renderer->AddActor(surfActor);
	renderer->SetViewport(0.0,0.0,0.5,1);///
	renderer->SetBackground(1.0, 1.0, 0.5);

////////////////////////////////show picture2//////////////////////////////////
	typedef itk::ImageToVTKImageFilter< FloatImageType> itkTovtkFilterType;
	itkTovtkFilterType::Pointer itkTovtkImageFilter2 = itkTovtkFilterType::New();
	itkTovtkImageFilter2->SetInput(DWI_castFilter->GetOutput());
	itkTovtkImageFilter2->Update();

	vtkSmartPointer<vtkMarchingCubes> ref_surface =
	    vtkSmartPointer<vtkMarchingCubes>::New();
	ref_surface->SetInputData(itkTovtkImageFilter2->GetOutput());
	ref_surface->ComputeNormalsOn();
	ref_surface->SetValue(0, isoValue);
		vtkSmartPointer<vtkPolyDataMapper> ref_surfMapper =
	    vtkSmartPointer<vtkPolyDataMapper>::New();
	ref_surfMapper->SetInputConnection(ref_surface->GetOutputPort());
	ref_surfMapper->ScalarVisibilityOff();

	vtkSmartPointer<vtkActor> ref_surfActor =
	    vtkSmartPointer<vtkActor>::New();
	ref_surfActor->SetMapper(ref_surfMapper);
	ref_surfActor->GetProperty()->SetColor(0.0, 1.0, 0.0);

	vtkSmartPointer<vtkRenderer> ref_renderer =
	    vtkSmartPointer<vtkRenderer>::New();
	ref_renderer->AddActor(ref_surfActor);
	ref_renderer->SetViewport(0.5,0.0,1,1);
	ref_renderer->SetBackground(1.0, 0.0, 1.0);

	vtkSmartPointer<vtkRenderWindow> renderWindow =
	    vtkSmartPointer<vtkRenderWindow>::New();
	renderWindow->AddRenderer(renderer);
	renderWindow->AddRenderer(ref_renderer);
	renderWindow->SetSize(640, 320);
	//renderWindow->Render();
	renderWindow->SetWindowName("局部配准选取初始点");
	renderWindow->Render();
	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
	    vtkSmartPointer<vtkRenderWindowInteractor>::New();//提供了一种独立于品台的交互举止相应鼠标案件始终等消息，监听到消息会自动调用回调函数。通过调用invokeevent函数将平台实践翻译成vtk事件

	renderWindowInteractor->SetRenderWindow(renderWindow);

	

	vtkSmartPointer<vtkRenderer> ref_renderer1 =
	    vtkSmartPointer<vtkRenderer>::New();
	ref_renderer1->AddActor(ref_surfActor);
	ref_renderer1->SetViewport(0.0,0.0,1,1);
	ref_renderer1->SetBackground(1.0, 0.0, 1.0);

	typedef itk::ResampleImageFilter<
	FloatImageType,
	FloatImageType >    ResampleFilterType;
	
	TransformType::Pointer finalTransform = TransformType::New();
	/*ifstream fin;
	string txtptn = "C:\\Users\\Administrator\\Desktop\\result\\zhangmengmeng\\1st\\2.txt";
	fin.open(txtptn.data());
	
	string s;
	while (getline(fin,s))
	{
		if(s=="Offset1")
		{
			getline(fin,s);//nextline
			
		}
	}*/

//	finalTransform->SetMatrix(,);
	//finalTransform->SetOffset(,,);
	ResampleFilterType::Pointer resampler = ResampleFilterType::New();
	resampler->SetTransform( finalTransform );
	resampler->SetInput( T1_castFilter->GetOutput() );/////////////全局图像


	vtkSmartPointer<vtkRenderWindow> renderWindowresult1 =
	    vtkSmartPointer<vtkRenderWindow>::New();
	renderWindowresult1->AddRenderer(ref_renderer1);//fixpicture
	renderWindowresult1->SetSize(640, 640);
	renderWindowresult1->SetWindowName("髁突配准的效果");
	renderWindowresult1->Render();
	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor1 =
	    vtkSmartPointer<vtkRenderWindowInteractor>::New();//提供了一种独立于品台的交互举止相应鼠标案件始终等消息，监听到消息会自动调用回调函数。通过调用invokeevent函数将平台实践翻译成vtk事件

	renderWindowInteractor1->SetRenderWindow(renderWindowresult1);
	//renderWindowInteractor1->Start();
	//renderWindowInteractor->Start();

	return EXIT_SUCCESS;
}
