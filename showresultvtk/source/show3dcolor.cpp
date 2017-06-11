/**********************************************************************

  �ļ���: 6.7_PolyDataMarchingCubes.cpp
  Copyright (c) ������, �޻���. All rights reserved.
  ������Ϣ�����: 
    http://www.vtkchina.org (VTK�й�)
	http://blog.csdn.net/www_doling_net (���鹤����) 

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
//���ԣ�../data/HeadMRVolume.mhd 200
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
	ReaderType::Pointer T1reader = ReaderType::New();//����һ��T1reader��ָ�����ڶ�ͼ
	typedef itk::GDCMImageIO       ImageIOType;
	ImageIOType::Pointer T1dicomIO = ImageIOType::New();///����T1dicomIOָ��
	T1reader->SetImageIO( T1dicomIO );//ʵ��ִ�ж�ȡ�����ImageIO�������ڱ����ӵ�ImageSeriesReader������ȷ������������������Ҫ��ȡ���ļ����͵�ImageIO������ȫ����
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
		//SeriesIdContainer���͵�T1seriesUID�洢������Ϣ����Ψһ��ʶDICOM��׼�и��ֲ�ͬ��Ϣ����
		const SeriesIdContainer & T1seriesUID = T1nameGenerator->GetSeriesUIDs();
		SeriesIdContainer::const_iterator T1seriesItr = T1seriesUID.begin();////������
		SeriesIdContainer::const_iterator T1seriesEnd = T1seriesUID.end();
		//��ʾ��ȡ����UID
		while( T1seriesItr != T1seriesEnd )
		{
			std::cout << T1seriesItr->c_str() << std::endl;
			++T1seriesItr;
		}
		std::string T1seriesIdentifier;

		T1seriesIdentifier = T1seriesUID.begin()->c_str();//ͨ����������ȡ���е�����Ƭ
		std::cout << std::endl << std::endl;
		std::cout << "Now reading series: " << std::endl << std::endl;
		std::cout << T1seriesIdentifier << std::endl;
		std::cout << std::endl << std::endl;
		FileNamesContainer T1fileNames;
		T1fileNames = T1nameGenerator->GetFileNames( T1seriesIdentifier );///��ȡ�ļ�
		T1reader->SetFileNames( T1fileNames );
		//���������writer��Update( )�����ܵ���ִ�С���ʱͼ�����Ƭ�����ڵ������ļ��У�ÿ���ļ�����һ����������Ƭ��������Щ��Ƭ���ļ������ļ�������������
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

	////////////////////////////////ʹ��itk���صڶ�������//////////////////////
	ReaderType::Pointer T2reader = ReaderType::New();
	ImageIOType::Pointer T2dicomIO = ImageIOType::New();
	T2reader->SetImageIO( T2dicomIO );
	NamesGeneratorType::Pointer T2nameGenerator = NamesGeneratorType::New();//��ȡʱ���ļ������ļ�������������
	T2nameGenerator->SetUseSeriesDetails( true );
//	T2nameGenerator->AddSeriesRestriction("0008|0021" );
	/*	T2nameGenerator->SetDirectory( "D:\\anonymized_before\\0020\\3" );*/
	//T2nameGenerator->SetDirectory( "C:\\Users\\Administrator\\Desktop\\devide\\haolinlin500\\condyle\\12haolinlinSEGandFILL"/*"C:\\Users\\Administrator\\Desktop\\devide\\square2" */);
	T2nameGenerator->SetDirectory(confixpth);
	//T2nameGenerator->SetDirectory( "C:\\Users\\Administrator\\Desktop\\niuxiahe\\mha�ļ���ȡ\\�½��ļ���\\2" );
	try
	{
		std::cout << std::endl << "The directory: " << std::endl;
		std::cout << std::endl << "16chengSEG"<< std::endl << std::endl;
		std::cout << "Contains the following DICOM Series: ";///ÿ��dicom���е���Ϣ��UID�ṩ
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
	typedef itk::Image< FloatPixelType, Dimension >      FloatImageType;////�����Ѿ���������������
	typedef itk::CastImageFilter<
	ImageType, FloatImageType >  CastFilterType;
//CastFilterType::Pointer       allT2_castFilter       = CastFilterType::New();
CastFilterType::Pointer       DWI_castFilter       = CastFilterType::New();
CastFilterType::Pointer       T1_castFilter       = CastFilterType::New();
		T1_castFilter->SetInput(       T1reader->GetOutput() );///����ȡ��������T1readerͨ��castfiltertype�õ��˲����ͼ��
	T1_castFilter->Update();

	//CastFilterType::Pointer       DWI_castFilter       = CastFilterType::New();
	DWI_castFilter->SetInput(       T2reader->GetOutput() );///��read��ͼ��Ϊ����
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
	DWI_filter->SetOutputMaximum( 255 );///�鵽0-255
	DWI_filter->SetInput( DWI_castFilter->GetOutput());
	DWI_filter->Update();
	typedef itk::ImageToVTKImageFilter< FloatImageType> itkTovtkFilterType;
	itkTovtkFilterType::Pointer itkTovtkImageFilter1 = itkTovtkFilterType::New();
	itkTovtkImageFilter1->SetInput(T1_castFilter->GetOutput());
	itkTovtkImageFilter1->Update();

	double isoValue = 255;

	vtkSmartPointer<vtkMarchingCubes> surface =
	    vtkSmartPointer<vtkMarchingCubes>::New();//�����ݽ���MarchingCubes �㷨�Ĵ�������ؽ�������������ݡ�
	surface->SetInputData(itkTovtkImageFilter1->GetOutput());
	surface->ComputeNormalsOn();//�����ֵ�淨���������Ⱦ����
	surface->SetValue(0, isoValue); 
  
	vtkSmartPointer<vtkPolyDataMapper> surfMapper =
	    vtkSmartPointer<vtkPolyDataMapper>::New();
	surfMapper->SetInputConnection(surface->GetOutputPort());//��Ⱦ����μ�������
	surfMapper->ScalarVisibilityOff();

	vtkSmartPointer<vtkActor> surfActor =
	    vtkSmartPointer<vtkActor>::New();//��Ⱦ���������ݵĿ��ӻ������ʾͼ��
	surfActor->SetMapper(surfMapper);
	surfActor->GetProperty()->SetColor(0.0, 0.0, 1.0);

	vtkSmartPointer<vtkRenderer> renderer =
	   vtkSmartPointer<vtkRenderer>::New();//�������������Ⱦ���̽����������������
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
	renderWindow->SetWindowName("�ֲ���׼ѡȡ��ʼ��");
	renderWindow->Render();
	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
	    vtkSmartPointer<vtkRenderWindowInteractor>::New();//�ṩ��һ�ֶ�����Ʒ̨�Ľ�����ֹ��Ӧ��갸��ʼ�յ���Ϣ����������Ϣ���Զ����ûص�������ͨ������invokeevent������ƽ̨ʵ�������vtk�¼�

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
	resampler->SetInput( T1_castFilter->GetOutput() );/////////////ȫ��ͼ��


	vtkSmartPointer<vtkRenderWindow> renderWindowresult1 =
	    vtkSmartPointer<vtkRenderWindow>::New();
	renderWindowresult1->AddRenderer(ref_renderer1);//fixpicture
	renderWindowresult1->SetSize(640, 640);
	renderWindowresult1->SetWindowName("��ͻ��׼��Ч��");
	renderWindowresult1->Render();
	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor1 =
	    vtkSmartPointer<vtkRenderWindowInteractor>::New();//�ṩ��һ�ֶ�����Ʒ̨�Ľ�����ֹ��Ӧ��갸��ʼ�յ���Ϣ����������Ϣ���Զ����ûص�������ͨ������invokeevent������ƽ̨ʵ�������vtk�¼�

	renderWindowInteractor1->SetRenderWindow(renderWindowresult1);
	//renderWindowInteractor1->Start();
	//renderWindowInteractor->Start();

	return EXIT_SUCCESS;
}
