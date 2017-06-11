/*#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#ifdef __BORLANDC__
#define ITK_LEAN_AND_MEAN
#endif*/

//#include "ITK2VTKP.h"



#include "itkGDCMImageIO.h"
#include "itkGDCMSeriesFileNames.h"
#include "itkImageSeriesReader.h"
#include "itkImageSeriesWriter.h"
// Software Guide : EndCodeSnippet

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
#include "vtkImageAnisotropicDiffusion3D.h"
#include "vtkImageCast.h"
#include <itkShiftScaleImageFilter.h>
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkPiecewiseFunction.h"
#include "vtkColorTransferFunction.h"
#include "vtkVolumeProperty.h"
#include "vtkVolumeRayCastIsosurfaceFunction.h"
#include "vtkVolumeRayCastCompositeFunction.h"
#include "vtkVolumeRayCastMapper.h"
#include "vtkVolume.h"
#include "vtkImageCast.h"
#include <vtkInteractorStyleImage.h>
#include <vtkDataSetMapper.h>
#include "vtkStripper.h"  
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



#include <vtkRendererCollection.h>
#include <vtkPointPicker.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkObjectFactory.h>



#include <vtkPoints.h>
#include <vtkLandmarkTransform.h>
#include <vtkMatrix4x4.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkVertexGlyphFilter.h>
#include <vtkAxesActor.h>


#include <vtkCoordinate.h>

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
ofstream fin;
string txtpth = "C:\\Users\\Administrator\\Desktop\\result\\2.txt";
string confixpth = "C:\\Users\\Administrator\\Desktop\\devide\\zhangmengmeng500\\concon\\8zhangmengmeng500SEGandFILL";
string conmovepth = "C:\\Users\\Administrator\\Desktop\\devide\\zhangmengmeng500\\concon\\2zhangnmengmeng500SEGandFILL";
string ludifixpth = "C:\\Users\\Administrator\\Desktop\\devide\\zhangmengmeng500\\ludi\\8zhangmengmeng";
string ludimovepth = "C:\\Users\\Administrator\\Desktop\\devide\\zhangmengmeng500\\ludi\\2zhangmengmeng";
string allfixpth = "C:\\Users\\Administrator\\Desktop\\niuxiahe\\maruohan-data\\Zhang Mengmeng (R)\\Zhang Mengmeng (R)\\zhangmengmeng20160817R\\ECIX3BKE\\CC1ZTGXR";
string allmovepth = "C:\\Users\\Administrator\\Desktop\\niuxiahe\\maruohan-data\\Zhang Mengmeng (R)\\Zhang Mengmeng (R)\\zhangmengmeng20160226R\\OU5R03VY\\WHPEAETD";
const char * savedicompth1 = "C:\\Users\\Administrator\\Desktop\\result\\condiff";
const char * savedicompth2 = "C:\\Users\\Administrator\\Desktop\\result\\conres";
const char * savedicompth3 = "C:\\Users\\Administrator\\Desktop\\result\\ludidiff";
const char * savedicompth4 = "C:\\Users\\Administrator\\Desktop\\result\\ludires";
const char * savedicompth5 = "C:\\Users\\Administrator\\Desktop\\result\\oridiff";

class CommandIterationUpdate : public itk::Command
{
public:
	typedef  CommandIterationUpdate   Self;
	typedef  itk::Command             Superclass;
	typedef itk::SmartPointer<Self>   Pointer;
	itkNewMacro( Self );

protected:
	CommandIterationUpdate() {};

public:
	typedef itk::RegularStepGradientDescentOptimizerv4<double>  OptimizerType;
	typedef   const OptimizerType *                             OptimizerPointer;
	void Execute(itk::Object *caller, const itk::EventObject & event) ITK_OVERRIDE
	{
		Execute( (const itk::Object *)caller, event);
	}
	void Execute(const itk::Object * object, const itk::EventObject & event) ITK_OVERRIDE
	{
		OptimizerPointer optimizer = static_cast< OptimizerPointer >( object );
		if( ! itk::IterationEvent().CheckEvent( &event ) )
		{
			return;
		}
		std::cout << optimizer->GetCurrentIteration() << "   ";
		std::cout << optimizer->GetValue() << "   ";
		std::cout << optimizer->GetCurrentPosition() << std::endl;

		
		fin.open(txtpth,ios::app);
		if(fin != 0)
		{
			fin<<"ѭ������"<<  optimizer->GetCurrentIteration()<<endl;
			fin<<  optimizer->GetValue() <<endl;
			fin<<  optimizer->GetCurrentPosition() <<endl;
			fin.close();
		}
	}//////////////���ÿ�ε�����Ĳ��ͼ��
};


	vtkSmartPointer<vtkRenderer> renderer =
		vtkSmartPointer<vtkRenderer>::New();

	vtkSmartPointer<vtkRenderer> ref_renderer =
	    vtkSmartPointer<vtkRenderer>::New();
	vtkSmartPointer<vtkRenderWindow> renderWindow =
		vtkSmartPointer<vtkRenderWindow>::New();

	vtkSmartPointer<vtkRenderer> rendererludi =
		vtkSmartPointer<vtkRenderer>::New();

	vtkSmartPointer<vtkRenderer> ref_rendererludi =
	    vtkSmartPointer<vtkRenderer>::New();
	vtkSmartPointer<vtkRenderWindow> renderWindowludi =
		vtkSmartPointer<vtkRenderWindow>::New();

	const unsigned int      Dimension = 3;
	typedef  float            FloatPixelType;
	typedef itk::Image< FloatPixelType, Dimension >  FixedImageType;
	typedef itk::Image< FloatPixelType, Dimension >  MovingImageType;

	typedef itk::Image< FloatPixelType, Dimension >      FloatImageType;////�����Ѿ���������������
	/*typedef itk::ThresholdImageFilter< FloatImageType >  ThresholdFilterType;  ///����˲�������ʹ�����ֲ�ͬ�ķ�ʽ��ת��һ��ͼ������ȼ�
	ThresholdFilterType::Pointer T1_filter = ThresholdFilterType::New();	
	ThresholdFilterType::Pointer DWI_filter = ThresholdFilterType::New();*/


	typedef itk::RescaleIntensityImageFilter<
	FloatImageType, FloatImageType > FloatRescaleFilterType;


FloatRescaleFilterType::Pointer DWI_filter = FloatRescaleFilterType::New();
	FloatRescaleFilterType::Pointer T1_filter = FloatRescaleFilterType::New();
	int num = 0;//num = 0.0;
	double x;
	double y;
	double z;
	double M00;
	double M01;
	double M02;
	double M10;
	double M11;
	double M12;
	double M20;
	double M21;
	double M22;
	double axis0,axis1,axis2;
	double rowangle;
	float a,b,c,w;
	string pth = "CC:\\Users\\Administrator\\Desktop\\samefiletext\\";
	

	typedef signed short    PixelType;
  typedef itk::Image< PixelType, Dimension >      ImageType;
 typedef itk::ImageSeriesReader< ImageType >     ReaderType;
 ReaderType::Pointer readerf = ReaderType::New();
 ReaderType::Pointer reader = ReaderType::New();
  typedef itk::ResampleImageFilter<
	  ImageType,
	ImageType >    ResampleFilterTypedicom;

	ResampleFilterTypedicom::Pointer resamplerdicom = ResampleFilterTypedicom::New();//�任���������׼ͼ��ӳ�䵽����ͼ��ռ�
	
class PointPickerInteractorStyle : public vtkInteractorStyleTrackballCamera
{
public:
	static double picked[8][3];//ÿ��ͼ��ѡȡ�ĸ���õ�������������
	static int pick_counter;

	static PointPickerInteractorStyle* New();
	vtkTypeMacro(PointPickerInteractorStyle, vtkInteractorStyleTrackballCamera);

	virtual void OnRightButtonDown()
	{
		std::cout << "Picking pixel: " << this->Interactor->GetEventPosition()[0] << " " << this->Interactor->GetEventPosition()[1] << std::endl;//���������Ļ����������Ϊ��λ��2��Ϊz����ͨ��Ϊ�㡾3��Ϊ����
		//vtkRendererͬ��Ҳ�����������view���꣨�����ͼ����Ⱦ����ϵͳ����display coordinates���豸����ʵ��screen���֮꣩��ִ������任
		vtkCollectionSimpleIterator temp;
		vtkRenderer *render;
		this->Interactor->GetRenderWindow()->GetRenderers()->InitTraversal(temp);///��������������������
		int *winSizePtr = this->Interactor->GetRenderWindow()->GetSize();///????????????�õ����ڴ�С
		if(pick_counter%2==0)///��������˫���¼��ȴ�������ͼ��ʼѡ��ÿ��ѡ�ĸ���
		{
			this->Interactor->GetPicker()->Pick(this->Interactor->GetEventPosition()[0],
			                                    this->Interactor->GetEventPosition()[1],
			                                    0,  // always zero.
			                                    this->Interactor->GetRenderWindow()->GetRenderers()->GetNextRenderer(temp));///��������vtkrender����
		}
		else//������
		{
 			vtkRendererCollection *rendercollection;
			
  			rendercollection=this->Interactor->GetRenderWindow()->GetRenderers();
  			rendercollection->GetNextRenderer(temp);
			render=rendercollection->GetNextRenderer(temp);
			this->Interactor->GetPicker()->Pick(this->Interactor->GetEventPosition()[0],
			                                    this->Interactor->GetEventPosition()[1],
			                                    0,  // always zero.
												render);
			                                   /* this->Interactor->GetRenderWindow()->GetRenderers()->GetNextRenderer(temp));*/

		}
		this->Interactor->GetPicker()->GetPickPosition(picked[pick_counter]);//�õ���������λ��pickedpositionΪ��������

		std::cout << "Picked value: " << picked[pick_counter][0] << " " << picked[pick_counter][1] << " " << picked[pick_counter][2] << std::endl;

		vtkSmartPointer<vtkSphereSource> sphereSource =
		    vtkSmartPointer<vtkSphereSource>::New();//ͨ���㷨ֱ����������
		sphereSource->Update();

		vtkSmartPointer<vtkPolyDataMapper> mapper =
		    vtkSmartPointer<vtkPolyDataMapper>::New();
		mapper->SetInputConnection(sphereSource->GetOutputPort());
		vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
		actor->SetMapper(mapper);


	//����ת��


		if(pick_counter%2==0)
		{
			actor->SetPosition(picked[pick_counter]);
		}
		else
		{
/*			picked[pick_counter][0]=picked[pick_counter][0]-val[0];*/
			actor->SetPosition(picked[pick_counter]);
		}
		actor->SetScale(0.5);
		actor->GetProperty()->SetColor(1.0, 0.0, 0.0);
		if(pick_counter%2==0)
		{
			this->Interactor->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->AddActor(actor);
		}
		else
		{
			render->AddActor(actor);
			/*this->Interactor->GetRenderWindow()->GetRenderers()->GetNextRenderer(temp)->AddActor(actor);*/
		}

		/*vtkInteractorStyleTrackballCamera::OnLeftButtonDown();*/
		pick_counter++;//û��һ�μ�һ
	}

	virtual void OnMiddleButtonDown()
	{
		vtkSmartPointer<vtkPoints> sourcePoints =
		    vtkSmartPointer<vtkPoints>::New();
		//double sourcePoint1[3] = {0.5, 0.0, 0.0};
		sourcePoints->InsertNextPoint(picked[0]);
		//double sourcePoint2[3] = {0.0, 0.5, 0.0};
		sourcePoints->InsertNextPoint(picked[2]);
		//double sourcePoint3[3] = {0.0, 0.0, 0.5};
		sourcePoints->InsertNextPoint(picked[4]);
		sourcePoints->InsertNextPoint(picked[6]);//���Ϊ����ͼ��


		vtkSmartPointer<vtkPoints> targetPoints =
		    vtkSmartPointer<vtkPoints>::New();
		//double targetPoint1[3] = {0.0, 0.0, 0.55};
		targetPoints->InsertNextPoint(picked[1]);
		//double targetPoint2[3] = {0.0, 0.55, 0.0};
		targetPoints->InsertNextPoint(picked[3]);
		//double targetPoint3[3] = {-0.55, 0.0, 0.0};
		targetPoints->InsertNextPoint(picked[5]);
		targetPoints->InsertNextPoint(picked[7]);///�Ҳ�ΪĿ��ͼ��
		
		vtkSmartPointer<vtkLandmarkTransform> landmarkTransform =
		    vtkSmartPointer<vtkLandmarkTransform>::New();//��ǵ���׼�������缫����׼���Ʒ��������С
		landmarkTransform->SetSourceLandmarks(sourcePoints);
		landmarkTransform->SetTargetLandmarks(targetPoints);
		landmarkTransform->SetModeToRigidBody();//������׼����Ϊ����///////////////////////���Ը�Ϊ�����任����
		landmarkTransform->Update();

		vtkSmartPointer<vtkPolyData> source =
		    vtkSmartPointer<vtkPolyData>::New();
		source->SetPoints(sourcePoints);

		vtkSmartPointer<vtkPolyData> target =
		    vtkSmartPointer<vtkPolyData>::New();
		target->SetPoints(targetPoints);

		vtkSmartPointer<vtkVertexGlyphFilter> sourceGlyphFilter =
		    vtkSmartPointer<vtkVertexGlyphFilter>::New();//��ʾ����׼�ĵ�
		sourceGlyphFilter->SetInputData(source);
		sourceGlyphFilter->Update();

		vtkSmartPointer<vtkVertexGlyphFilter> targetGlyphFilter =
		    vtkSmartPointer<vtkVertexGlyphFilter>::New();//��ʾĿ���
		targetGlyphFilter->SetInputData(target);
		targetGlyphFilter->Update();

		vtkSmartPointer<vtkTransformPolyDataFilter> transformFilter =
		    vtkSmartPointer<vtkTransformPolyDataFilter>::New();//��Ԫ��ǵ���б任����ʾ��׼��ĵ㼯
		transformFilter->SetInputData(sourceGlyphFilter->GetOutput());//��Դ��任����ʾ
		transformFilter->SetTransform(landmarkTransform);
		transformFilter->Update();
		//std::cout << "landmarkTransform �� " << landmarkTransform->GetReferenceCount()<<  std::endl;


		vtkSmartPointer<vtkPolyDataMapper> solutionMapper =
		vtkSmartPointer<vtkPolyDataMapper>::New();
	solutionMapper->SetInputConnection(transformFilter->GetOutputPort());

	vtkSmartPointer<vtkActor> solutionActor =
		vtkSmartPointer<vtkActor>::New();
	solutionActor->SetMapper(solutionMapper);
	solutionActor->GetProperty()->SetColor(0.2,0.5,1);/////////////////�ı����￴���ܹ�����ɫ�任
	solutionActor->GetProperty()->SetPointSize(9);

	ref_renderer->AddActor(solutionActor);
	renderWindow->AddRenderer(ref_renderer);
	renderWindow->Render();

	landmarkTransform->Inverse();
	vtkSmartPointer<vtkMatrix4x4> minv = landmarkTransform->GetMatrix() ; 
	//minv = landmarkTransform->GetElements(0,0);
	std::cout << "VTK�õ��ĳ�ʼ����The resulting inverse matrix is: " << *(minv) << std::endl; 
	  /*for(int i = 0;i<= 3;i++)
    {
		std::cout <<" "<<std::endl;
        for(int j = 0;j <= 3;j++)
        {
           std::cout<<minv->Element[i][j]<<std::endl;
			
        }
	  }*/
  ////�õ���ʼ�任λ�õ���itk�е���׼����registration������ͼ�񣬱任������
 ///�õ�����///
	   x = minv->Element[0][3];
	   y = minv->Element[1][3];
	   z = minv->Element[2][3];
	   M00 = minv->Element[0][0];
	   M01 = minv->Element[0][1];
	   M02 = minv->Element[0][2];
	   M10 = minv->Element[1][0];
	   M11 = minv->Element[1][1];
	   M12 = minv->Element[1][2];
	   M20 = minv->Element[2][0];
	   M21 = minv->Element[2][1];
	   M22 = minv->Element[2][2];

	  std::cout<<"VTK��ʼƽ�ƣ�"<<x<<" "<<y<<" "<<z<<" "<<std::endl;

	  ori_transform();
	  fin.open(txtpth,ios::app);
		if(fin != 0)
		{
			fin << "��ʼ��ת"<<axis0<<","<<axis1<<","<<axis2<<","<<rowangle<<endl;
			fin<<"��ʼƽ�ƣ�"<<x<<" "<<y<<" "<<z<<" "<<endl;
			fin.close();
		}
	 // allregistion();
	 registion();
	  allregistion();
	}

};

vtkStandardNewMacro(PointPickerInteractorStyle);
int PointPickerInteractorStyle::pick_counter = 0;
double PointPickerInteractorStyle::picked[8][3]= {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
typedef itk::CastImageFilter<
	ImageType, FloatImageType >  CastFilterType;
CastFilterType::Pointer       allT2_castFilter       = CastFilterType::New();
CastFilterType::Pointer       DWI_castFilter       = CastFilterType::New();
CastFilterType::Pointer       T1_castFilter       = CastFilterType::New();
int main()
{
	
	/*typedef signed short    Pixel;	
	typedef itk::Image< Pixel, Dimension >         ImageType;
	typedef itk::ImageSeriesReader< ImageType >        ReaderType;*/
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


	//seg();
	////////////ͼ��Ԥ����
	/*typedef   float            FloatPixelType;
	typedef itk::Image< FloatPixelType, Dimension >   FloatImageType;////�����Ѿ���������������*/
	///ʹ�ö����ͼ�����������˲�������ʵ����
	//typedef itk::CastImageFilter<
	//ImageType, FloatImageType >  CastFilterType;
	//ͨ������New( )���������������˲����������ָ��itk::SmartPointers

	//CastFilterType::Pointer       T1_castFilter       = CastFilterType::New();
	T1_castFilter->SetInput(       T1reader->GetOutput() );///����ȡ��������T1readerͨ��castfiltertype�õ��˲����ͼ��
	T1_castFilter->Update();

	//CastFilterType::Pointer       DWI_castFilter       = CastFilterType::New();
	DWI_castFilter->SetInput(       T2reader->GetOutput() );///��read��ͼ��Ϊ����
	DWI_castFilter->Update();

	/*typedef itk::RescaleIntensityImageFilter<
	FloatImageType, FloatImageType > FloatRescaleFilterType;



	FloatRescaleFilterType::Pointer T1_rescalercast = FloatRescaleFilterType::New();*/

	T1_filter->SetOutputMinimum(   0 );
	T1_filter->SetOutputMaximum( 255 );
	T1_filter->SetInput( T1_castFilter->GetOutput());
	T1_filter->Update();
	//FloatRescaleFilterType::Pointer DWI_rescalercast = FloatRescaleFilterType::New();
	DWI_filter->SetOutputMinimum(   0 );
	DWI_filter->SetOutputMaximum( 255 );///�鵽0-255
	DWI_filter->SetInput( DWI_castFilter->GetOutput());
	DWI_filter->Update();

	//typedef itk::ThresholdImageFilter< FloatImageType >  ThresholdFilterType;  ///����˲�������ʹ�����ֲ�ͬ�ķ�ʽ��ת��һ��ͼ������ȼ�
	//ThresholdFilterType::Pointer T1_filter = ThresholdFilterType::New();
	/*T1_filter->SetInput( T1_rescalercast->GetOutput() );
	T1_filter->ThresholdOutside( 0,255 );////������ѡ��
	T1_filter->Update();////T1filter�еõ�Ҫ��ʾ��ͼ��

	//ThresholdFilterType::Pointer DWI_filter = ThresholdFilterType::New();
	DWI_filter->SetInput( DWI_rescalercast->GetOutput() );
	DWI_filter->ThresholdOutside( 0,255 );
	DWI_filter->Update();*/

	///���任���ӵ���׼������


	//ITK TO VTK
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

	//vtkSmartPointer<vtkRenderer> renderer =
	  //  vtkSmartPointer<vtkRenderer>::New();//�������������Ⱦ���̽����������������
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

	//vtkSmartPointer<vtkRenderer> ref_renderer =
	  //  vtkSmartPointer<vtkRenderer>::New();
	ref_renderer->AddActor(ref_surfActor);
	ref_renderer->SetViewport(0.5,0.0,1,1);
	ref_renderer->SetBackground(1.0, 0.0, 1.0);

	//vtkSmartPointer<vtkRenderWindow> renderWindow =
	  //  vtkSmartPointer<vtkRenderWindow>::New();
	renderWindow->AddRenderer(renderer);
	renderWindow->AddRenderer(ref_renderer);
	renderWindow->SetSize(640, 320);
	//renderWindow->Render();
	renderWindow->SetWindowName("�ֲ���׼ѡȡ��ʼ��");
	renderWindow->Render();
	vtkSmartPointer<vtkPointPicker> pointPicker =
	    vtkSmartPointer<vtkPointPicker>::New();

	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
	    vtkSmartPointer<vtkRenderWindowInteractor>::New();//�ṩ��һ�ֶ�����Ʒ̨�Ľ�����ֹ��Ӧ��갸��ʼ�յ���Ϣ����������Ϣ���Զ����ûص�������ͨ������invokeevent������ƽ̨ʵ�������vtk�¼�
	renderWindowInteractor->SetPicker(pointPicker);
	renderWindowInteractor->SetRenderWindow(renderWindow);

	vtkSmartPointer<PointPickerInteractorStyle> style =
	    vtkSmartPointer<PointPickerInteractorStyle>::New();
	renderWindowInteractor->SetInteractorStyle( style );
	renderWindowInteractor->Start();
	
	return EXIT_SUCCESS;

}

void ori_transform()//��ת��������Ԫ����ת��
{
	float fourw = M00+M11+M22;
	float foura = M00-M11-M22;
	float fourb = M11-M00-M22;
	float fourc = M22-M00-M11;

	///�Ƚϴ�С
	int big = 0;//���w���
	float fourf = fourw;
	if(foura>fourf)
	{
		fourf = foura;
		big = 1;
	}
	if(fourb>fourf)
	{
		fourf = fourb;
		big = 2;
	}
	if(fourc>fourf)
	{
		fourf = fourc;
		big = 3;
	}

	float bigV = sqrt(fourf+1.0f)*0.5f;
	float mult = 0.25f/bigV;

	//����big�ж����յļ��㹫ʽ
	switch(big)
	{
	case 0:
		w = bigV;
		a = (M12-M21)*mult;
		b = (M20-M02)*mult;
		c = (M01-M10)*mult;
		break;
	case 1:
		a = bigV;
		w = (M12-M21)*mult;
		b = (M10+M01)*mult;
		c = (M20+M02)*mult;
		break;
	case 2:
		b = bigV;
		w = (M20-M02)*mult;
		a = (M01+M10)*mult;
		c = (M12+M21)*mult;
		break;
	case 3:
		c = bigV;
		w = (M01-M10)*mult;
		a = (M20+M02)*mult;
		b = (M12+M21)*mult;
		break;
	}
	///��Ԫ��W,A,B,C���潫��Ԫ��ת��Ϊ��ת��ͽǶ�
	rowangle = acos(w)*2;
	double sina = sin(rowangle/2);
	axis0 = a/sina;
	axis1 = b/sina;
	axis2 = c/sina;

}
typedef itk::VersorRigid3DTransform< double > TransformType;
TransformType::Pointer finalTransform = TransformType::New();
typedef itk::ResampleImageFilter<
	MovingImageType,
	FixedImageType >    ResampleFilterType;

	ResampleFilterType::Pointer resamplermha = ResampleFilterType::New();//�任���������׼ͼ��ӳ�䵽����ͼ��ռ�
void registion()//�ֲ���׼
{

	/*double max;
	max = x;*/
	//typedef itk::VersorRigid3DTransform< double > TransformType;
	typedef itk::RegularStepGradientDescentOptimizerv4<double>    OptimizerType;//�ݶ��½�����һ��
	typedef itk::MeanSquaresImageToImageMetricv4<
	FixedImageType,
	MovingImageType >   MetricType;
	typedef itk::ImageRegistrationMethodv4<
	FixedImageType,
	MovingImageType,
	TransformType >           RegistrationType;

	MetricType::Pointer         metric        = MetricType::New();
	OptimizerType::Pointer      optimizer     = OptimizerType::New();
	RegistrationType::Pointer   registration  = RegistrationType::New();
	registration->SetMetric(        metric        );
	registration->SetOptimizer(     optimizer     );
	TransformType::Pointer  initialTransform = TransformType::New();
	
	typedef itk::ImageFileReader< FixedImageType  > FixedImageReaderType;
	typedef itk::ImageFileReader< MovingImageType > MovingImageReaderType;

	registration->SetFixedImage(    DWI_filter->GetOutput()    );
	registration->SetMovingImage(   T1_filter->GetOutput()   );

	typedef itk::CenteredTransformInitializer<
	TransformType,
	FixedImageType,
	MovingImageType >  TransformInitializerType;
	TransformInitializerType::Pointer initializer =
	    TransformInitializerType::New();
		initializer->SetTransform(   initialTransform );

	initializer->SetFixedImage(  DWI_filter->GetOutput() );
	initializer->SetMovingImage( T1_filter->GetOutput() );

	initializer->GeometryOn();//����MomentsOn( )��ѡ����������
	initializer->InitializeTransform();/////���������任���ĺ�ƽ�Ƶļ���

	///���ݳ�ʼ�任������Ԫ�����ڶ�����תq = w+xi+yj+zk///
	typedef TransformType::VersorType  VersorType;
	typedef VersorType::VectorType     VectorType;
	VersorType     rotation;
	VectorType     axis;//��ת��
	VectorType		trans;//ƽ����
	//if(num==0)//��һ�εĲ�����VTK����
	
	axis[0] = axis0;
	axis[1] = axis1;
	axis[2] = axis2;
	const double angle = rowangle;//���ݱ任����õ���ת�Ƕȣ�����������������������������������������������������
	rotation.Set(  axis, angle  );
	std::cout << "��VTK����õ���ʼ��ת��"<<axis0<<","<<axis1<<","<<axis2<<std::endl;
	initialTransform->SetRotation( rotation );
	trans[0] = x;
	trans[1] = y;
	trans[2] = z;
	
	
	initialTransform->SetTranslation(trans);
	////��γ�ʼ��ƽ����Ϣ
	registration->SetInitialTransform( initialTransform );
	 //�����Ż��������.  
    //ע��: ��ת�ͱ任��"��λ�̶�" �ǲ�ͬ��,��ת�û��ȶ���, ƽ���Ժ��׶��� 
	typedef OptimizerType::ScalesType       OptimizerScalesType;
	OptimizerScalesType optimizerScales( initialTransform->GetNumberOfParameters() );
	const double translationScale = 1.0/1000.0 ;//����λ��ת��ƽ�Ƶı�������ܴ����Ǿ����Ż����ṩ�ı��������Ե��ŵ����
	optimizerScales[0] = 1.0;
	optimizerScales[1] = 1.0;
	optimizerScales[2] = 1.0;
	optimizerScales[3] = translationScale;
	optimizerScales[4] = translationScale;
	optimizerScales[5] = translationScale;
	optimizer->SetScales( optimizerScales ); 


	optimizer->SetNumberOfIterations(100);///�����
	optimizer->SetLearningRate( 0.5 );//����ϵ��ѧϰ��?
	optimizer->SetMinimumStepLength( 0.001 );///������Ĺ���
	optimizer->SetReturnBestParametersAndValue(true);
	//optimizer->set

	CommandIterationUpdate::Pointer observer = CommandIterationUpdate::New();///��ת���������.ʵ����һ Command/Observer ����, ������׼���̵�ִ��, ��������׼���̵�ִ��.  
	optimizer->AddObserver( itk::IterationEvent(), observer );
	//////////////////��֪�����Ǹ�ʲô��
	const unsigned int numberOfLevels = 1;
	RegistrationType::ShrinkFactorsArrayType shrinkFactorsPerLevel;
	shrinkFactorsPerLevel.SetSize( 1 );
	shrinkFactorsPerLevel[0] = 1;

	RegistrationType::SmoothingSigmasArrayType smoothingSigmasPerLevel;
	smoothingSigmasPerLevel.SetSize( 1 );
	smoothingSigmasPerLevel[0] = 0;

	registration->SetNumberOfLevels( numberOfLevels );//���������ڽ������ж�ֱ��ʼ��Ĳ���
	registration->SetSmoothingSigmasPerLevel( smoothingSigmasPerLevel );
	registration->SetShrinkFactorsPerLevel( shrinkFactorsPerLevel );
	/////////////////////////////
	try
	{
		registration->Update();// //������׼���̵�ִ��
		std::cout << "Optimizer stop condition: "
		          << registration->GetOptimizer()->GetStopConditionDescription()
		          << std::endl;
		fin.open(txtpth,ios::app);
		if(fin != 0)
		{
			fin<<   "ֹͣ����"<<registration->GetOptimizer()->GetStopConditionDescription()<<endl;
			
			fin.close();
		}
	}
	catch( itk::ExceptionObject & err )
	{
		std::cerr << "ExceptionObject caught !" << std::endl;
		std::cerr << err << std::endl;
		
	}
	//��ȡ���ղ���
	const TransformType::ParametersType finalParameters =
	    registration->GetOutput()->Get()->GetParameters();

	const double versorX              = finalParameters[0];
	const double versorY              = finalParameters[1];
	const double versorZ              = finalParameters[2];
	const double finalTranslationX    = finalParameters[3];
	const double finalTranslationY    = finalParameters[4];
	const double finalTranslationZ    = finalParameters[5];
	const unsigned int numberOfIterations = optimizer->GetCurrentIteration();
	const double bestValue = optimizer->GetValue();

	// Print out results
	//
	std::cout << std::endl << std::endl;
	std::cout << "���ݵ�ITK��ı任Result = " << std::endl;
	std::cout << " versor X      = " << versorX  << std::endl;
	std::cout << " versor Y      = " << versorY  << std::endl;
	std::cout << " versor Z      = " << versorZ  << std::endl;
	std::cout << " Translation X = " << finalTranslationX  << std::endl;
	std::cout << " Translation Y = " << finalTranslationY  << std::endl;
	std::cout << " Translation Z = " << finalTranslationZ  << std::endl;
	std::cout << " Iterations    = " << numberOfIterations << std::endl;
	std::cout << " Metric value  = " << bestValue          << std::endl;

	
	//TransformType::Pointer finalTransform = TransformType::New();
	//��6��������ȥ��ʵ�ʵ���ת�����ƫ��������˵����
	const TransformType::ParametersType finalParametersfix =
			registration->GetOutput()->Get()->GetFixedParameters() ;
	finalTransform->SetFixedParameters( registration->GetOutput()->Get()->GetFixedParameters() );
	finalTransform->SetParameters( finalParameters );

	// Software Guide : BeginCodeSnippet
	TransformType::MatrixType matrix = finalTransform->GetMatrix();
	TransformType::OffsetType offset = finalTransform->GetOffset();
	
	std::cout << "ITK��תMatrix = " << std::endl << matrix << std::endl;//��ת����
	std::cout << "ITKƽ��Offset = " << std::endl << offset << std::endl;//ƫ����
	fin.open(txtpth,ios::app);
		if(fin != 0)
		{
			fin<<"finalParametersfix"<<finalParametersfix<<endl;
			
			fin << "��ֹ��תMatrix = " << std::endl << matrix << endl;
			fin<< "��ֹƽ��Offset = " << std::endl << offset << endl;
			fin.close();
		}
	//showresult();
	//�ֲ���������ȫ��
	x = offset[0];
	y = offset[1];
	z = offset[2];
	M00 = matrix[0][0];
	M01 =matrix[0][1];
	M02 = matrix[0][2];
	M10 = matrix[1][0];
	M11 = matrix[1][1];
	M12 = matrix[1][2];
	M20 =matrix[2][0];
	M21 = matrix[2][1];
	M22 =matrix[2][2];
	ori_transform();
	std::cout << "���ݵ�����ı仯����"<<axis0<<","<<axis1<<","<<axis2<<","<<rowangle<<std::endl;
	fin.open(txtpth,ios::app);
		if(fin != 0)
		{
			fin << "��ֹ��ת"<<axis0<<","<<axis1<<","<<axis2<<","<<rowangle<<endl;
			fin.close();
		}
		fin.open(txtpth,ios::app);
		if(fin != 0)
		{
			std::cout << std::endl << std::endl;
			fin << "���ݵ�ITK��ı任Result = " << std::endl;
			fin<< " versor X      = " << versorX  << std::endl;
			fin << " versor Y      = " << versorY  << std::endl;
			fin << " versor Z      = " << versorZ  << std::endl;
			fin << " Translation X = " << finalTranslationX  << std::endl;
			fin << " Translation Y = " << finalTranslationY  << std::endl;
			fin<< " Translation Z = " << finalTranslationZ  << std::endl;
			fin << " Iterations    = " << numberOfIterations << std::endl;
			fin << " Metric value  = " << bestValue          << std::endl;
			fin.close();
		}
	/*typedef itk::ResampleImageFilter<
	MovingImageType,
	FixedImageType >    ResampleFilterType;
	saveresulttodicom();
	ResampleFilterType::Pointer resamplermha = ResampleFilterType::New();//�任���������׼ͼ��ӳ�䵽����ͼ��ռ�*/
	
	saveresulttomha();
	saveresulttodicom();
	if(num==1)
	{
		showresult();
	}
	num++;
}


vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor1 =
	    vtkSmartPointer<vtkRenderWindowInteractor>::New();//�ṩ��һ�ֶ�����Ʒ̨�Ľ�����ֹ��Ӧ��갸��ʼ�յ���Ϣ����������Ϣ���Զ����ûص�������ͨ������invokeevent������ƽ̨ʵ�������vtk�¼�
void allregistion()//ȫ����׼���벢����ȫ��ͼ��
{
	//////////////////�����һ��ȫ��ͼ��///////////////////////////
	typedef signed short    Pixel;
	
	typedef itk::Image< Pixel, Dimension >         ImageType;
	typedef itk::ImageSeriesReader< ImageType >        ReaderType;
	ReaderType::Pointer allT1reader = ReaderType::New();//����һ��T1reader��ָ�����ڶ�ͼ
	typedef itk::GDCMImageIO       ImageIOType;
	ImageIOType::Pointer allT1dicomIO = ImageIOType::New();///����T1dicomIOָ��
	allT1reader->SetImageIO( allT1dicomIO );//ʵ��ִ�ж�ȡ�����ImageIO�������ڱ����ӵ�ImageSeriesReader������ȷ������������������Ҫ��ȡ���ļ����͵�ImageIO������ȫ����
	typedef itk::GDCMSeriesFileNames NamesGeneratorType;
	NamesGeneratorType::Pointer allT1nameGenerator = NamesGeneratorType::New();
	allT1nameGenerator->SetUseSeriesDetails( true );
	allT1nameGenerator->SetDirectory(ludimovepth);
	
	typedef std::vector< std::string >   FileNamesContainer;
	try
	{
		std::cout << std::endl << "The directory: " << std::endl;
		std::cout << std::endl << "first"<< std::endl << std::endl;
		std::cout << "Contains the following DICOM Series: ";
		std::cout << std::endl << std::endl;
		typedef std::vector< std::string >    SeriesIdContainer;
		const SeriesIdContainer & allT1seriesUID = allT1nameGenerator->GetSeriesUIDs();
		SeriesIdContainer::const_iterator allT1seriesItr = allT1seriesUID.begin();////������
		SeriesIdContainer::const_iterator allT1seriesEnd = allT1seriesUID.end();
		//��ʾ��ȡ����UID
		while( allT1seriesItr != allT1seriesEnd )
		{
			std::cout << allT1seriesItr->c_str() << std::endl;
			++allT1seriesItr;
		}
		std::string allT1seriesIdentifier;

		allT1seriesIdentifier = allT1seriesUID.begin()->c_str();//ͨ����������ȡ���е�����Ƭ
		std::cout << std::endl << std::endl;
		std::cout << "Now reading series: " << std::endl << std::endl;
		std::cout << allT1seriesIdentifier << std::endl;
		std::cout << std::endl << std::endl;
		FileNamesContainer allT1fileNames;
		allT1fileNames = allT1nameGenerator->GetFileNames( allT1seriesIdentifier );///��ȡ�ļ�
		allT1reader->SetFileNames( allT1fileNames );
		try
		{
			allT1reader->Update();
		}
		catch (itk::ExceptionObject &ex)
		{
			std::cout << ex << std::endl;
			//return 1;
		}

	}
	catch (itk::ExceptionObject &ex)
	{
		std::cout << ex << std::endl;
		//return EXIT_FAILURE;
	}
	/////////////////����ڶ���//////////////////
		typedef signed short    Pixel;
	
	typedef itk::Image< Pixel, Dimension >         ImageType;
	typedef itk::ImageSeriesReader< ImageType >        ReaderType;
	ReaderType::Pointer allT2reader = ReaderType::New();//����һ��T1reader��ָ�����ڶ�ͼ
	typedef itk::GDCMImageIO       ImageIOType;
	ImageIOType::Pointer allT2dicomIO = ImageIOType::New();///����T1dicomIOָ��
	allT2reader->SetImageIO( allT2dicomIO );//ʵ��ִ�ж�ȡ�����ImageIO�������ڱ����ӵ�ImageSeriesReader������ȷ������������������Ҫ��ȡ���ļ����͵�ImageIO������ȫ����
	typedef itk::GDCMSeriesFileNames NamesGeneratorType;
	NamesGeneratorType::Pointer allT2nameGenerator = NamesGeneratorType::New();
	allT2nameGenerator->SetUseSeriesDetails( true );
	allT2nameGenerator->SetDirectory(ludifixpth);
	
	typedef std::vector< std::string >   FileNamesContainer;
	try
	{
		std::cout << std::endl << "The directory: " << std::endl;
		std::cout << std::endl << "second:"<< std::endl << std::endl;
		std::cout << "Contains the following DICOM Series: ";
		std::cout << std::endl << std::endl;
		typedef std::vector< std::string >    SeriesIdContainer;
		const SeriesIdContainer & allT2seriesUID = allT2nameGenerator->GetSeriesUIDs();
		SeriesIdContainer::const_iterator allT2seriesItr = allT2seriesUID.begin();////������
		SeriesIdContainer::const_iterator allT2seriesEnd = allT2seriesUID.end();
		//��ʾ��ȡ����UID
		while( allT2seriesItr != allT2seriesEnd )
		{
			std::cout << allT2seriesItr->c_str() << std::endl;
			++allT2seriesItr;
		}
		std::string allT2seriesIdentifier;

		allT2seriesIdentifier = allT2seriesUID.begin()->c_str();//ͨ����������ȡ���е�����Ƭ
		std::cout << std::endl << std::endl;
		std::cout << "Now reading series: " << std::endl << std::endl;
		std::cout << allT2seriesIdentifier << std::endl;
		std::cout << std::endl << std::endl;
		FileNamesContainer allT2fileNames;
		allT2fileNames = allT2nameGenerator->GetFileNames( allT2seriesIdentifier );///��ȡ�ļ�
		allT2reader->SetFileNames( allT2fileNames );
		try
		{
			allT2reader->Update();
		}
		catch (itk::ExceptionObject &ex)
		{
			std::cout << ex << std::endl;
			//return 1;
		}

	}
	catch (itk::ExceptionObject &ex)
	{
		std::cout << ex << std::endl;
		//return EXIT_FAILURE;
	}
	////ͼ��Ԥ����
	//typedef itk::CastImageFilter<
	//ImageType, FloatImageType >  CastFilterType;
	//ͨ������New( )���������������˲����������ָ��itk::SmartPointers

	CastFilterType::Pointer       allT1_castFilter       = CastFilterType::New();
	allT1_castFilter->SetInput(       allT1reader->GetOutput() );///����ȡ��������T1readerͨ��castfiltertype�õ��˲����ͼ��
	allT1_castFilter->Update();


	
	allT2_castFilter->SetInput(       allT2reader->GetOutput() );///��read��ͼ��Ϊ����
	allT2_castFilter->Update();

	//typedef itk::RescaleIntensityImageFilter<
	//FloatImageType, FloatImageType > FloatRescaleFilterType;



	//FloatRescaleFilterType::Pointer T1_filter = FloatRescaleFilterType::New();

	T1_filter->SetOutputMinimum(   0 );
	T1_filter->SetOutputMaximum( 255 );
	T1_filter->SetInput( allT1_castFilter->GetOutput());
	T1_filter->Update();
	//FloatRescaleFilterType::Pointer DWI_filter = FloatRescaleFilterType::New();
	DWI_filter->SetOutputMinimum(   0 );
	DWI_filter->SetOutputMaximum( 255 );///�鵽0-255
	DWI_filter->SetInput( allT2_castFilter->GetOutput());
	
	DWI_filter->Update();


	//ITK TO VTK
	typedef itk::ImageToVTKImageFilter< FloatImageType> itkTovtkFilterType;
	itkTovtkFilterType::Pointer itkTovtkImageFilter1 = itkTovtkFilterType::New();
	itkTovtkImageFilter1->SetInput(allT1_castFilter->GetOutput());
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

	//vtkSmartPointer<vtkRenderer> renderer =
	  //  vtkSmartPointer<vtkRenderer>::New();//�������������Ⱦ���̽����������������
	rendererludi->AddActor(surfActor);
	rendererludi->SetViewport(0.0,0.0,0.5,1);///
	rendererludi->SetBackground(1.0, 1.0, 0.5);

////////////////////////////////show picture2//////////////////////////////////
	typedef itk::ImageToVTKImageFilter< FloatImageType> itkTovtkFilterType;
	itkTovtkFilterType::Pointer itkTovtkImageFilter2 = itkTovtkFilterType::New();
	itkTovtkImageFilter2->SetInput(allT2_castFilter->GetOutput());
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

	//vtkSmartPointer<vtkRenderer> ref_renderer =
	  //  vtkSmartPointer<vtkRenderer>::New();
	ref_rendererludi->AddActor(ref_surfActor);
	ref_rendererludi->SetViewport(0.5,0.0,1,1);
	ref_rendererludi->SetBackground(1.0, 0.0, 1.0);

	//vtkSmartPointer<vtkRenderWindow> renderWindow =
	  //  vtkSmartPointer<vtkRenderWindow>::New();
	renderWindowludi->AddRenderer(rendererludi);
	renderWindowludi->AddRenderer(ref_rendererludi);
	renderWindowludi->SetSize(640, 320);
	//renderWindowludi->Render();
	renderWindowludi->SetWindowName("­����׼");
	renderWindowludi->Render();
	
	renderWindowInteractor1->SetRenderWindow(renderWindowludi);

	vtkSmartPointer<PointPickerInteractorStyle> style =
	    vtkSmartPointer<PointPickerInteractorStyle>::New();
	renderWindowInteractor1->SetInteractorStyle( style );
	//renderWindowInteractor->Start();
	registion();
	/*vtkSmartPointer<vtkPointPicker> pointPicker =
	    vtkSmartPointer<vtkPointPicker>::New();

	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
	    vtkSmartPointer<vtkRenderWindowInteractor>::New();//�ṩ��һ�ֶ�����Ʒ̨�Ľ�����ֹ��Ӧ��갸��ʼ�յ���Ϣ����������Ϣ���Զ����ûص�������ͨ������invokeevent������ƽ̨ʵ�������vtk�¼�
	renderWindowInteractor->SetPicker(pointPicker);
	renderWindowInteractor->SetRenderWindow(renderWindow);

	vtkSmartPointer<PointPickerInteractorStyle> style =
	    vtkSmartPointer<PointPickerInteractorStyle>::New();
	renderWindowInteractor->SetInteractorStyle( style );
	renderWindowInteractor->Start();*/
	
	
}

void saveresulttodicom()//��׼�󽫽�������dicom
{
	////////////////////fix
  typedef itk::GDCMImageIO                        ImageIOType;
  typedef itk::GDCMSeriesFileNames                NamesGeneratorType;

  ImageIOType::Pointer gdcmIOf = ImageIOType::New();
  NamesGeneratorType::Pointer namesGeneratorf = NamesGeneratorType::New();
   //namesGeneratorf->SetInputDirectory("C:\\Users\\Administrator\\Desktop\\niuxiahe\\maruohan-data\\haoliangliang(left)\\haoliangliang20151211L\\ZC255IGE\\C0LUTGR3"/*"C:\\Users\\Administrator\\Desktop\\devide\\haolinlin500\\testsquare12"*/ );
  namesGeneratorf->SetInputDirectory(allfixpth);//ע�������Ƿ���
  const ReaderType::FileNamesContainer & filenamesf =
                            namesGeneratorf->GetInputFileNames();
  // Software Guide : EndCodeSnippet

  std::size_t numberOfFileNamesf = filenamesf.size();
  std::cout << numberOfFileNamesf << std::endl;

  //ReaderType::Pointer readerf = ReaderType::New();

  readerf->SetImageIO( gdcmIOf );
  readerf->SetFileNames( filenamesf );
  try
    {
    // Software Guide : BeginCodeSnippet
    readerf->Update();
    // Software Guide : EndCodeSnippet
    }
  catch (itk::ExceptionObject &excp)
    {
    std::cerr << "Exception throiting the image" << std::endl;
    std::cerr << excp << std::endl;
    //return EXIT_FAILURE;
    }

	///////////////////move
	//typedef signed short    PixelType;
 // const unsigned int      Dimension = 3;

 // typedef itk::Image< PixelType, Dimension >      ImageType;
  //typedef itk::ImageSeriesReader< ImageType >     ReaderType;
  typedef itk::GDCMImageIO                        ImageIOType;
  typedef itk::GDCMSeriesFileNames                NamesGeneratorType;

  ImageIOType::Pointer gdcmIO = ImageIOType::New();
  NamesGeneratorType::Pointer namesGenerator = NamesGeneratorType::New();
  namesGenerator->SetInputDirectory(allmovepth);

  const ReaderType::FileNamesContainer & filenames =
                            namesGenerator->GetInputFileNames();
  // Software Guide : EndCodeSnippet

  std::size_t numberOfFileNames = filenames.size();
  std::cout << numberOfFileNames << std::endl;

  reader->SetImageIO( gdcmIO );
  reader->SetFileNames( filenames );
  try
    {
    // Software Guide : BeginCodeSnippet
    reader->Update();
    // Software Guide : EndCodeSnippet
    }
  catch (itk::ExceptionObject &excp)
    {
    std::cerr << "Exception te" << std::endl;
    std::cerr << excp << std::endl;
    //return EXIT_FAILURE;
    }
 /* typedef itk::ResampleImageFilter<
	  ImageType,
	ImageType >    ResampleFilterType;

	ResampleFilterType::Pointer resampler = ResampleFilterType::New();//�任���������׼ͼ��ӳ�䵽����ͼ��ռ�*/

	resamplerdicom->SetTransform( finalTransform );
	resamplerdicom->SetInput( reader->GetOutput() );/////////////ȫ��ͼ��

	ImageType::Pointer fixedImagenew = readerf->GetOutput();

	resamplerdicom->SetSize(    fixedImagenew->GetLargestPossibleRegion().GetSize() );
	resamplerdicom->SetOutputOrigin(  fixedImagenew->GetOrigin() );
	resamplerdicom->SetOutputSpacing( fixedImagenew->GetSpacing() );
	resamplerdicom->SetOutputDirection( fixedImagenew->GetDirection() );
	  const char * outputDirectory;
	  if(num ==0)
	  {
		  outputDirectory = savedicompth2;
	  }
	  if(num!=0)
	  {
		  outputDirectory =savedicompth4;
	  }
  itksys::SystemTools::MakeDirectory( outputDirectory );
  typedef signed short    OutputPixelType;
  const unsigned int      OutputDimension = 2;

  typedef itk::Image< OutputPixelType, OutputDimension >    Image2DType;

  typedef itk::ImageSeriesWriter<
                             ImageType, Image2DType >  SeriesWriterType;
  SeriesWriterType::Pointer seriesWriter = SeriesWriterType::New();

  seriesWriter->SetInput( resamplerdicom->GetOutput() );
  seriesWriter->SetImageIO( gdcmIO );///////////////////////////   
  namesGenerator->SetOutputDirectory( outputDirectory );

  seriesWriter->SetFileNames( namesGenerator->GetOutputFileNames() );
  seriesWriter->SetMetaDataDictionaryArray(
                        reader->GetMetaDataDictionaryArray() );
   try
    {
    seriesWriter->Update();
    }
  catch( itk::ExceptionObject & excp )
    {
    std::cerr << "Exception thrown while writing the series " << std::endl;
    std::cerr << excp << std::endl;
    //return EXIT_FAILURE;
    }


  typedef itk::SubtractImageFilter<
	ImageType,
	ImageType,
	ImageType > DifferenceFilterType;//����
	DifferenceFilterType::Pointer difference = DifferenceFilterType::New();

	typedef itk::RescaleIntensityImageFilter<
	ImageType,
	ImageType >   RescalerType;///���µ���������ʹ���Ǹ�������ʹ��ֵ���ӻ���Ϊ����
	RescalerType::Pointer intensityRescaler = RescalerType::New();

	intensityRescaler->SetInput( difference->GetOutput() );
	intensityRescaler->SetOutputMinimum(   0 );
	intensityRescaler->SetOutputMaximum( 255 );

	difference->SetInput1( readerf->GetOutput() );//�̶�ͼ��
	difference->SetInput2( resamplerdicom->GetOutput() );//���²�����ĸ���ͼ��

	resamplerdicom->SetDefaultPixelValue( 1 );
	const char * outputDirectory2;
	if(num==0)
	{
		outputDirectory2 = savedicompth1;
	}
	if(num!=0)
	{
		outputDirectory2 =savedicompth3;
	}
  itksys::SystemTools::MakeDirectory( outputDirectory2);
  typedef signed short    OutputPixelType;
  //const unsigned int      OutputDimension = 2;

  typedef itk::Image< OutputPixelType, OutputDimension >    Image2DType;

  typedef itk::ImageSeriesWriter<
                             ImageType, Image2DType >  SeriesWriterType;
  SeriesWriterType::Pointer seriesWriter2 = SeriesWriterType::New();

  seriesWriter2->SetInput( difference->GetOutput() );
  seriesWriter2->SetImageIO( gdcmIOf );
  namesGeneratorf->SetOutputDirectory( outputDirectory2 );

  seriesWriter2->SetFileNames( namesGeneratorf->GetOutputFileNames() );
  seriesWriter2->SetMetaDataDictionaryArray(
                        readerf->GetMetaDataDictionaryArray() );
   try
    {
    seriesWriter2->Update();
    }
  catch( itk::ExceptionObject & excp )
    {
    std::cerr << "Exception thrown while writing the series " << std::endl;
    std::cerr << excp << std::endl;
    //return EXIT_FAILURE;
    }



  if(num==0)
  {

	DifferenceFilterType::Pointer differenceori = DifferenceFilterType::New();

	typedef itk::RescaleIntensityImageFilter<
	ImageType,
	ImageType >   RescalerType;///���µ���������ʹ���Ǹ�������ʹ��ֵ���ӻ���Ϊ����
	RescalerType::Pointer intensityRescalerori = RescalerType::New();

	intensityRescalerori->SetInput( differenceori->GetOutput() );
	intensityRescalerori->SetOutputMinimum(   0 );
	intensityRescalerori->SetOutputMaximum( 255 );

	differenceori->SetInput1( readerf->GetOutput() );//�̶�ͼ��
	differenceori->SetInput2( reader->GetOutput() );//���²�����ĸ���ͼ��

	resamplerdicom->SetDefaultPixelValue( 1 );
	const char * outputDirectory3;

	outputDirectory3 = savedicompth5;
	
  itksys::SystemTools::MakeDirectory( outputDirectory3);
  typedef signed short    OutputPixelType;
  //const unsigned int      OutputDimension = 2;

  typedef itk::Image< OutputPixelType, OutputDimension >    Image2DType;

  typedef itk::ImageSeriesWriter<
                             ImageType, Image2DType >  SeriesWriterType;
  SeriesWriterType::Pointer seriesWriter3 = SeriesWriterType::New();

  seriesWriter3->SetInput( differenceori->GetOutput() );
  seriesWriter3->SetImageIO( gdcmIOf );
  namesGeneratorf->SetOutputDirectory( outputDirectory3 );

  seriesWriter3->SetFileNames( namesGeneratorf->GetOutputFileNames() );
  seriesWriter3->SetMetaDataDictionaryArray(
                        readerf->GetMetaDataDictionaryArray() );
   try
    {
    seriesWriter3->Update();
    }
  catch( itk::ExceptionObject & excp )
    {
    std::cerr << "Exception thrown while writing the series " << std::endl;
    std::cerr << excp << std::endl;
    //return EXIT_FAILURE;
    }
  }
  //addcolor();
}


/*void addcolor()
{
		

	typedef itk::CastImageFilter<
	ImageType, FloatImageType >  CastFilterType;
	CastFilterType::Pointer       move_castFilter       = CastFilterType::New();
	move_castFilter->SetInput(       reader->GetOutput() );///����ȡ��������T1readerͨ��castfiltertype�õ��˲����ͼ��
	move_castFilter->Update();

	CastFilterType::Pointer       fix_castFilter       = CastFilterType::New();
	fix_castFilter->SetInput(       readerf->GetOutput() );///��read��ͼ��Ϊ����
	fix_castFilter->Update();

	typedef itk::ImageToVTKImageFilter< FloatImageType> itkTovtkFilterType;
	itkTovtkFilterType::Pointer moveimgvtk = itkTovtkFilterType::New();
	moveimgvtk->SetInput( move_castFilter->GetOutput() );
	moveimgvtk->Update();

	typedef itk::ImageToVTKImageFilter< FloatImageType> itkTovtkFilterType;
	itkTovtkFilterType::Pointer fiximgvtk = itkTovtkFilterType::New();
	fiximgvtk->SetInput(fix_castFilter->GetOutput());
	fiximgvtk->Update();
	//////////////////////////

	
}*/
vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractorresult =
	    vtkSmartPointer<vtkRenderWindowInteractor>::New();//�ṩ��һ�ֶ�����Ʒ̨�Ľ�����ֹ��Ӧ��갸��ʼ�յ���Ϣ����������Ϣ���Զ����ûص�������ͨ������invokeevent������ƽ̨ʵ�������vtk�¼�
void showresult()
{
		double isoValue = 255;
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

	vtkSmartPointer<vtkActor> ref_surfActorfix =
	    vtkSmartPointer<vtkActor>::New();
	ref_surfActorfix->SetMapper(ref_surfMapper);
	ref_surfActorfix->GetProperty()->SetColor(0.0, 1.0, 0.0);	
	
	///////////�任���ͼ��
	ResampleFilterType::Pointer resampler = ResampleFilterType::New();
	resampler->SetTransform( finalTransform );
	resampler->SetInput( T1_castFilter->GetOutput() );/////////////ȫ��ͼ��

	FixedImageType::Pointer fixedImage = DWI_castFilter->GetOutput();

	resampler->SetSize(    fixedImage->GetLargestPossibleRegion().GetSize() );
	resampler->SetOutputOrigin(  fixedImage->GetOrigin() );
	resampler->SetOutputSpacing( fixedImage->GetSpacing() );
	resampler->SetOutputDirection( fixedImage->GetDirection() );
	resampler->SetDefaultPixelValue( 100 );


	typedef itk::ImageToVTKImageFilter< FloatImageType> itkTovtkFilterType;
	itkTovtkFilterType::Pointer moveimgvtk = itkTovtkFilterType::New();
	moveimgvtk->SetInput( resampler->GetOutput() );
	moveimgvtk->Update();
	//moveimgvtk->GetOutput()->GetDimensions();

	vtkSmartPointer<vtkRenderWindow> renderWindowresultcon =
		vtkSmartPointer<vtkRenderWindow>::New();

	vtkSmartPointer<vtkRenderer> rendererresultcon =
		vtkSmartPointer<vtkRenderer>::New();
	

		vtkSmartPointer<vtkMarchingCubes> ref_surface2 =
	    vtkSmartPointer<vtkMarchingCubes>::New();
	ref_surface2->SetInputData(moveimgvtk->GetOutput());
	ref_surface2->ComputeNormalsOn();
	ref_surface2->SetValue(0, isoValue);
		vtkSmartPointer<vtkPolyDataMapper> ref_surfMapper2 =
	    vtkSmartPointer<vtkPolyDataMapper>::New();
	ref_surfMapper2->SetInputConnection(ref_surface2->GetOutputPort());
	ref_surfMapper2->ScalarVisibilityOff();

	vtkSmartPointer<vtkActor> ref_surfActorres =
	    vtkSmartPointer<vtkActor>::New();
	ref_surfActorres->SetMapper(ref_surfMapper2);
	ref_surfActorres->GetProperty()->SetColor(0.0, 0.0, 1.0);

	
	rendererresultcon->AddActor(ref_surfActorfix);
	rendererresultcon->SetViewport(0.0,0.0,1,1);
	rendererresultcon->SetBackground(1.0, 0.0, 1.0);

	renderWindowresultcon->AddRenderer(rendererresultcon);
	//renderWindow->AddRenderer(ref_renderer);
	renderWindowresultcon->SetSize(640, 640);
	renderWindowresultcon->SetWindowName("�ֲ���׼Ч��");
	rendererresultcon->AddActor(ref_surfActorres);
	renderWindowresultcon->AddRenderer(rendererresultcon);
	renderWindowresultcon->Render();
	
	renderWindowInteractorresult->SetRenderWindow(renderWindowresultcon);

	vtkSmartPointer<PointPickerInteractorStyle> style =
	    vtkSmartPointer<PointPickerInteractorStyle>::New();
	renderWindowInteractorresult->SetInteractorStyle( style );
	
		//renderWindowInteractor->Start();
		//renderWindowInteractor1->Start();
	
}
void saveresulttomha()
{
		


	resamplermha->SetTransform( finalTransform );
	resamplermha->SetInput( T1_filter->GetOutput() );/////////////ȫ��ͼ��

	FixedImageType::Pointer fixedImage = DWI_filter->GetOutput();

	resamplermha->SetSize(    fixedImage->GetLargestPossibleRegion().GetSize() );
	resamplermha->SetOutputOrigin(  fixedImage->GetOrigin() );
	resamplermha->SetOutputSpacing( fixedImage->GetSpacing() );
	resamplermha->SetOutputDirection( fixedImage->GetDirection() );
	resamplermha->SetDefaultPixelValue( 100 );

	typedef  unsigned char                                          OutputPixelType;
	typedef itk::Image< OutputPixelType, Dimension >                OutputImageType;
	typedef itk::CastImageFilter< FixedImageType, OutputImageType > CastFilterTypeN;//ת���ز�����ͼ����������͵�����writer����������
	typedef itk::ImageFileWriter< OutputImageType >                 WriterType;

	WriterType::Pointer      writer =  WriterType::New();
	CastFilterTypeN::Pointer  caster =  CastFilterTypeN::New();
	string snum = to_string(num); 
	string cnum = pth+"deform" + snum+".mha";
	writer->SetFileName( cnum );

	caster->SetInput( resamplermha->GetOutput() );
	writer->SetInput( caster->GetOutput()   );
	//writer->Update();
}