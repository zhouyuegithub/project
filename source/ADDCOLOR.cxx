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
			fin<<"循环过程"<<  optimizer->GetCurrentIteration()<<endl;
			fin<<  optimizer->GetValue() <<endl;
			fin<<  optimizer->GetCurrentPosition() <<endl;
			fin.close();
		}
	}//////////////添加每次迭代后的差分图像
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

	typedef itk::Image< FloatPixelType, Dimension >      FloatImageType;////输入已经定义过，这是输出
	/*typedef itk::ThresholdImageFilter< FloatImageType >  ThresholdFilterType;  ///这个滤波器可以使用三种不同的方式来转化一个图像的亮度级
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

	ResampleFilterTypedicom::Pointer resamplerdicom = ResampleFilterTypedicom::New();//变换结果将待配准图像映射到参照图像空间
	
class PointPickerInteractorStyle : public vtkInteractorStyleTrackballCamera
{
public:
	static double picked[8][3];//每侧图像选取四个点得到的是世界坐标
	static int pick_counter;

	static PointPickerInteractorStyle* New();
	vtkTypeMacro(PointPickerInteractorStyle, vtkInteractorStyleTrackballCamera);

	virtual void OnRightButtonDown()
	{
		std::cout << "Picking pixel: " << this->Interactor->GetEventPosition()[0] << " " << this->Interactor->GetEventPosition()[1] << std::endl;//鼠标点击的屏幕坐标以像素为单位【2】为z坐标通常为零【3】为对象
		//vtkRenderer同样也在世界坐标和view坐标（计算机图形渲染坐标系统）和display coordinates（设备上事实的screen坐标）之间执行坐标变换
		vtkCollectionSimpleIterator temp;
		vtkRenderer *render;
		this->Interactor->GetRenderWindow()->GetRenderers()->InitTraversal(temp);///？？？？？？？？？？
		int *winSizePtr = this->Interactor->GetRenderWindow()->GetSize();///????????????得到窗口大小
		if(pick_counter%2==0)///如果点击了双数下即先从左侧的视图开始选点每侧选四个点
		{
			this->Interactor->GetPicker()->Pick(this->Interactor->GetEventPosition()[0],
			                                    this->Interactor->GetEventPosition()[1],
			                                    0,  // always zero.
			                                    this->Interactor->GetRenderWindow()->GetRenderers()->GetNextRenderer(temp));///窗口坐标vtkrender对象
		}
		else//单数下
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
		this->Interactor->GetPicker()->GetPickPosition(picked[pick_counter]);//得到世界坐标位置pickedposition为世界坐标

		std::cout << "Picked value: " << picked[pick_counter][0] << " " << picked[pick_counter][1] << " " << picked[pick_counter][2] << std::endl;

		vtkSmartPointer<vtkSphereSource> sphereSource =
		    vtkSmartPointer<vtkSphereSource>::New();//通过算法直接生成数据
		sphereSource->Update();

		vtkSmartPointer<vtkPolyDataMapper> mapper =
		    vtkSmartPointer<vtkPolyDataMapper>::New();
		mapper->SetInputConnection(sphereSource->GetOutputPort());
		vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
		actor->SetMapper(mapper);


	//坐标转换


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
		pick_counter++;//没点一次加一
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
		sourcePoints->InsertNextPoint(picked[6]);//左侧为浮动图像


		vtkSmartPointer<vtkPoints> targetPoints =
		    vtkSmartPointer<vtkPoints>::New();
		//double targetPoint1[3] = {0.0, 0.0, 0.55};
		targetPoints->InsertNextPoint(picked[1]);
		//double targetPoint2[3] = {0.0, 0.55, 0.0};
		targetPoints->InsertNextPoint(picked[3]);
		//double targetPoint3[3] = {-0.55, 0.0, 0.0};
		targetPoints->InsertNextPoint(picked[5]);
		targetPoints->InsertNextPoint(picked[7]);///右侧为目标图像
		
		vtkSmartPointer<vtkLandmarkTransform> landmarkTransform =
		    vtkSmartPointer<vtkLandmarkTransform>::New();//标记点配准法两个电极在配准后的品均距离最小
		landmarkTransform->SetSourceLandmarks(sourcePoints);
		landmarkTransform->SetTargetLandmarks(targetPoints);
		landmarkTransform->SetModeToRigidBody();//设置配准类型为刚体///////////////////////可以改为其他变换类型
		landmarkTransform->Update();

		vtkSmartPointer<vtkPolyData> source =
		    vtkSmartPointer<vtkPolyData>::New();
		source->SetPoints(sourcePoints);

		vtkSmartPointer<vtkPolyData> target =
		    vtkSmartPointer<vtkPolyData>::New();
		target->SetPoints(targetPoints);

		vtkSmartPointer<vtkVertexGlyphFilter> sourceGlyphFilter =
		    vtkSmartPointer<vtkVertexGlyphFilter>::New();//显示代配准的点
		sourceGlyphFilter->SetInputData(source);
		sourceGlyphFilter->Update();

		vtkSmartPointer<vtkVertexGlyphFilter> targetGlyphFilter =
		    vtkSmartPointer<vtkVertexGlyphFilter>::New();//显示目标点
		targetGlyphFilter->SetInputData(target);
		targetGlyphFilter->Update();

		vtkSmartPointer<vtkTransformPolyDataFilter> transformFilter =
		    vtkSmartPointer<vtkTransformPolyDataFilter>::New();//对元标记点进行变换来显示配准后的点集
		transformFilter->SetInputData(sourceGlyphFilter->GetOutput());//将源点变换后显示
		transformFilter->SetTransform(landmarkTransform);
		transformFilter->Update();
		//std::cout << "landmarkTransform ： " << landmarkTransform->GetReferenceCount()<<  std::endl;


		vtkSmartPointer<vtkPolyDataMapper> solutionMapper =
		vtkSmartPointer<vtkPolyDataMapper>::New();
	solutionMapper->SetInputConnection(transformFilter->GetOutputPort());

	vtkSmartPointer<vtkActor> solutionActor =
		vtkSmartPointer<vtkActor>::New();
	solutionActor->SetMapper(solutionMapper);
	solutionActor->GetProperty()->SetColor(0.2,0.5,1);/////////////////改变这里看看能够由颜色变换
	solutionActor->GetProperty()->SetPointSize(9);

	ref_renderer->AddActor(solutionActor);
	renderWindow->AddRenderer(ref_renderer);
	renderWindow->Render();

	landmarkTransform->Inverse();
	vtkSmartPointer<vtkMatrix4x4> minv = landmarkTransform->GetMatrix() ; 
	//minv = landmarkTransform->GetElements(0,0);
	std::cout << "VTK得到的初始矩阵The resulting inverse matrix is: " << *(minv) << std::endl; 
	  /*for(int i = 0;i<= 3;i++)
    {
		std::cout <<" "<<std::endl;
        for(int j = 0;j <= 3;j++)
        {
           std::cout<<minv->Element[i][j]<<std::endl;
			
        }
	  }*/
  ////得到初始变换位置调用itk中的配准部分registration（两幅图像，变换参数）
 ///得到参数///
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

	  std::cout<<"VTK初始平移："<<x<<" "<<y<<" "<<z<<" "<<std::endl;

	  ori_transform();
	  fin.open(txtpth,ios::app);
		if(fin != 0)
		{
			fin << "初始旋转"<<axis0<<","<<axis1<<","<<axis2<<","<<rowangle<<endl;
			fin<<"初始平移："<<x<<" "<<y<<" "<<z<<" "<<endl;
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


	//seg();
	////////////图像预处理
	/*typedef   float            FloatPixelType;
	typedef itk::Image< FloatPixelType, Dimension >   FloatImageType;////输入已经定义过，这是输出*/
	///使用定义的图像类型来对滤波器进行实例化
	//typedef itk::CastImageFilter<
	//ImageType, FloatImageType >  CastFilterType;
	//通过调用New( )操作来创建对象滤波器并将结果指向itk::SmartPointers

	//CastFilterType::Pointer       T1_castFilter       = CastFilterType::New();
	T1_castFilter->SetInput(       T1reader->GetOutput() );///将读取到的数据T1reader通过castfiltertype得到滤波后的图像
	T1_castFilter->Update();

	//CastFilterType::Pointer       DWI_castFilter       = CastFilterType::New();
	DWI_castFilter->SetInput(       T2reader->GetOutput() );///将read的图作为输入
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
	DWI_filter->SetOutputMaximum( 255 );///归到0-255
	DWI_filter->SetInput( DWI_castFilter->GetOutput());
	DWI_filter->Update();

	//typedef itk::ThresholdImageFilter< FloatImageType >  ThresholdFilterType;  ///这个滤波器可以使用三种不同的方式来转化一个图像的亮度级
	//ThresholdFilterType::Pointer T1_filter = ThresholdFilterType::New();
	/*T1_filter->SetInput( T1_rescalercast->GetOutput() );
	T1_filter->ThresholdOutside( 0,255 );////参数的选择
	T1_filter->Update();////T1filter中得到要显示的图像

	//ThresholdFilterType::Pointer DWI_filter = ThresholdFilterType::New();
	DWI_filter->SetInput( DWI_rescalercast->GetOutput() );
	DWI_filter->ThresholdOutside( 0,255 );
	DWI_filter->Update();*/

	///将变换连接到配准类型中


	//ITK TO VTK
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

	//vtkSmartPointer<vtkRenderer> renderer =
	  //  vtkSmartPointer<vtkRenderer>::New();//负责管理场景的渲染过程将场景对象组合起来
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
	renderWindow->SetWindowName("局部配准选取初始点");
	renderWindow->Render();
	vtkSmartPointer<vtkPointPicker> pointPicker =
	    vtkSmartPointer<vtkPointPicker>::New();

	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
	    vtkSmartPointer<vtkRenderWindowInteractor>::New();//提供了一种独立于品台的交互举止相应鼠标案件始终等消息，监听到消息会自动调用回调函数。通过调用invokeevent函数将平台实践翻译成vtk事件
	renderWindowInteractor->SetPicker(pointPicker);
	renderWindowInteractor->SetRenderWindow(renderWindow);

	vtkSmartPointer<PointPickerInteractorStyle> style =
	    vtkSmartPointer<PointPickerInteractorStyle>::New();
	renderWindowInteractor->SetInteractorStyle( style );
	renderWindowInteractor->Start();
	
	return EXIT_SUCCESS;

}

void ori_transform()//旋转矩阵与四元数的转化
{
	float fourw = M00+M11+M22;
	float foura = M00-M11-M22;
	float fourb = M11-M00-M22;
	float fourc = M22-M00-M11;

	///比较大小
	int big = 0;//如果w最大
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

	//根据big判断最终的计算公式
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
	///四元数W,A,B,C下面将四元数转化为旋转轴和角度
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

	ResampleFilterType::Pointer resamplermha = ResampleFilterType::New();//变换结果将待配准图像映射到参照图像空间
void registion()//局部配准
{

	/*double max;
	max = x;*/
	//typedef itk::VersorRigid3DTransform< double > TransformType;
	typedef itk::RegularStepGradientDescentOptimizerv4<double>    OptimizerType;//梯度下降法的一种
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

	initializer->GeometryOn();//访问MomentsOn( )来选择力矩中心
	initializer->InitializeTransform();/////用来触发变换中心和平移的计算

	///传递初始变换参数四元数用于定义旋转q = w+xi+yj+zk///
	typedef TransformType::VersorType  VersorType;
	typedef VersorType::VectorType     VectorType;
	VersorType     rotation;
	VectorType     axis;//旋转轴
	VectorType		trans;//平移量
	//if(num==0)//第一次的参数由VTK传递
	
	axis[0] = axis0;
	axis[1] = axis1;
	axis[2] = axis2;
	const double angle = rowangle;//根据变换矩阵得到旋转角度？？？？？？？？？？？？？？？？？？？？？？？？？？？
	rotation.Set(  axis, angle  );
	std::cout << "由VTK计算得到初始旋转轴"<<axis0<<","<<axis1<<","<<axis2<<std::endl;
	initialTransform->SetRotation( rotation );
	trans[0] = x;
	trans[1] = y;
	trans[2] = z;
	
	
	initialTransform->SetTranslation(trans);
	////如何初始化平移信息
	registration->SetInitialTransform( initialTransform );
	 //设置优化组件参数.  
    //注意: 旋转和变换的"单位刻度" 是不同的,旋转用弧度度量, 平移以毫米度量 
	typedef OptimizerType::ScalesType       OptimizerScalesType;
	OptimizerScalesType optimizerScales( initialTransform->GetNumberOfParameters() );
	const double translationScale = 1.0/1000.0 ;//：单位旋转和平移的比例差异很大，我们就以优化器提供的比例范函性的优点举例
	optimizerScales[0] = 1.0;
	optimizerScales[1] = 1.0;
	optimizerScales[2] = 1.0;
	optimizerScales[3] = translationScale;
	optimizerScales[4] = translationScale;
	optimizerScales[5] = translationScale;
	optimizer->SetScales( optimizerScales ); 


	optimizer->SetNumberOfIterations(100);///最大数
	optimizer->SetLearningRate( 0.5 );//阻尼系数学习率?
	optimizer->SetMinimumStepLength( 0.001 );///收敛域的公差
	optimizer->SetReturnBestParametersAndValue(true);
	//optimizer->set

	CommandIterationUpdate::Pointer observer = CommandIterationUpdate::New();///跳转后进行运算.实例化一 Command/Observer 对象, 监视配准过程的执行, 并触发配准过程的执行.  
	optimizer->AddObserver( itk::IterationEvent(), observer );
	//////////////////不知道这是干什么的
	const unsigned int numberOfLevels = 1;
	RegistrationType::ShrinkFactorsArrayType shrinkFactorsPerLevel;
	shrinkFactorsPerLevel.SetSize( 1 );
	shrinkFactorsPerLevel[0] = 1;

	RegistrationType::SmoothingSigmasArrayType smoothingSigmasPerLevel;
	smoothingSigmasPerLevel.SetSize( 1 );
	smoothingSigmasPerLevel[0] = 0;

	registration->SetNumberOfLevels( numberOfLevels );//用于设置在金字塔中多分辨率级的层数
	registration->SetSmoothingSigmasPerLevel( smoothingSigmasPerLevel );
	registration->SetShrinkFactorsPerLevel( shrinkFactorsPerLevel );
	/////////////////////////////
	try
	{
		registration->Update();// //触发配准过程的执行
		std::cout << "Optimizer stop condition: "
		          << registration->GetOptimizer()->GetStopConditionDescription()
		          << std::endl;
		fin.open(txtpth,ios::app);
		if(fin != 0)
		{
			fin<<   "停止条件"<<registration->GetOptimizer()->GetStopConditionDescription()<<endl;
			
			fin.close();
		}
	}
	catch( itk::ExceptionObject & err )
	{
		std::cerr << "ExceptionObject caught !" << std::endl;
		std::cerr << err << std::endl;
		
	}
	//获取最终参数
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
	std::cout << "传递到ITK后的变换Result = " << std::endl;
	std::cout << " versor X      = " << versorX  << std::endl;
	std::cout << " versor Y      = " << versorY  << std::endl;
	std::cout << " versor Z      = " << versorZ  << std::endl;
	std::cout << " Translation X = " << finalTranslationX  << std::endl;
	std::cout << " Translation Y = " << finalTranslationY  << std::endl;
	std::cout << " Translation Z = " << finalTranslationZ  << std::endl;
	std::cout << " Iterations    = " << numberOfIterations << std::endl;
	std::cout << " Metric value  = " << bestValue          << std::endl;

	
	//TransformType::Pointer finalTransform = TransformType::New();
	//从6个参数中去看实际的旋转矩阵和偏移量最有说服力
	const TransformType::ParametersType finalParametersfix =
			registration->GetOutput()->Get()->GetFixedParameters() ;
	finalTransform->SetFixedParameters( registration->GetOutput()->Get()->GetFixedParameters() );
	finalTransform->SetParameters( finalParameters );

	// Software Guide : BeginCodeSnippet
	TransformType::MatrixType matrix = finalTransform->GetMatrix();
	TransformType::OffsetType offset = finalTransform->GetOffset();
	
	std::cout << "ITK旋转Matrix = " << std::endl << matrix << std::endl;//旋转矩阵
	std::cout << "ITK平移Offset = " << std::endl << offset << std::endl;//偏移量
	fin.open(txtpth,ios::app);
		if(fin != 0)
		{
			fin<<"finalParametersfix"<<finalParametersfix<<endl;
			
			fin << "终止旋转Matrix = " << std::endl << matrix << endl;
			fin<< "终止平移Offset = " << std::endl << offset << endl;
			fin.close();
		}
	//showresult();
	//局部参数传递全局
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
	std::cout << "传递到整体的变化参数"<<axis0<<","<<axis1<<","<<axis2<<","<<rowangle<<std::endl;
	fin.open(txtpth,ios::app);
		if(fin != 0)
		{
			fin << "终止旋转"<<axis0<<","<<axis1<<","<<axis2<<","<<rowangle<<endl;
			fin.close();
		}
		fin.open(txtpth,ios::app);
		if(fin != 0)
		{
			std::cout << std::endl << std::endl;
			fin << "传递到ITK后的变换Result = " << std::endl;
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
	ResampleFilterType::Pointer resamplermha = ResampleFilterType::New();//变换结果将待配准图像映射到参照图像空间*/
	
	saveresulttomha();
	saveresulttodicom();
	if(num==1)
	{
		showresult();
	}
	num++;
}


vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor1 =
	    vtkSmartPointer<vtkRenderWindowInteractor>::New();//提供了一种独立于品台的交互举止相应鼠标案件始终等消息，监听到消息会自动调用回调函数。通过调用invokeevent函数将平台实践翻译成vtk事件
void allregistion()//全局配准读入并处理全局图像
{
	//////////////////读入第一幅全局图像///////////////////////////
	typedef signed short    Pixel;
	
	typedef itk::Image< Pixel, Dimension >         ImageType;
	typedef itk::ImageSeriesReader< ImageType >        ReaderType;
	ReaderType::Pointer allT1reader = ReaderType::New();//定义一个T1reader的指针用于读图
	typedef itk::GDCMImageIO       ImageIOType;
	ImageIOType::Pointer allT1dicomIO = ImageIOType::New();///定义T1dicomIO指针
	allT1reader->SetImageIO( allT1dicomIO );//实际执行读取任务的ImageIO对象现在被连接到ImageSeriesReader。这是确保我们用适于我们想要读取的文件类型的ImageIO对象的最安全方法
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
		SeriesIdContainer::const_iterator allT1seriesItr = allT1seriesUID.begin();////迭代器
		SeriesIdContainer::const_iterator allT1seriesEnd = allT1seriesUID.end();
		//显示读取到的UID
		while( allT1seriesItr != allT1seriesEnd )
		{
			std::cout << allT1seriesItr->c_str() << std::endl;
			++allT1seriesItr;
		}
		std::string allT1seriesIdentifier;

		allT1seriesIdentifier = allT1seriesUID.begin()->c_str();//通过迭代器读取所有单张切片
		std::cout << std::endl << std::endl;
		std::cout << "Now reading series: " << std::endl << std::endl;
		std::cout << allT1seriesIdentifier << std::endl;
		std::cout << std::endl << std::endl;
		FileNamesContainer allT1fileNames;
		allT1fileNames = allT1nameGenerator->GetFileNames( allT1seriesIdentifier );///读取文件
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
	/////////////////读入第二幅//////////////////
		typedef signed short    Pixel;
	
	typedef itk::Image< Pixel, Dimension >         ImageType;
	typedef itk::ImageSeriesReader< ImageType >        ReaderType;
	ReaderType::Pointer allT2reader = ReaderType::New();//定义一个T1reader的指针用于读图
	typedef itk::GDCMImageIO       ImageIOType;
	ImageIOType::Pointer allT2dicomIO = ImageIOType::New();///定义T1dicomIO指针
	allT2reader->SetImageIO( allT2dicomIO );//实际执行读取任务的ImageIO对象现在被连接到ImageSeriesReader。这是确保我们用适于我们想要读取的文件类型的ImageIO对象的最安全方法
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
		SeriesIdContainer::const_iterator allT2seriesItr = allT2seriesUID.begin();////迭代器
		SeriesIdContainer::const_iterator allT2seriesEnd = allT2seriesUID.end();
		//显示读取到的UID
		while( allT2seriesItr != allT2seriesEnd )
		{
			std::cout << allT2seriesItr->c_str() << std::endl;
			++allT2seriesItr;
		}
		std::string allT2seriesIdentifier;

		allT2seriesIdentifier = allT2seriesUID.begin()->c_str();//通过迭代器读取所有单张切片
		std::cout << std::endl << std::endl;
		std::cout << "Now reading series: " << std::endl << std::endl;
		std::cout << allT2seriesIdentifier << std::endl;
		std::cout << std::endl << std::endl;
		FileNamesContainer allT2fileNames;
		allT2fileNames = allT2nameGenerator->GetFileNames( allT2seriesIdentifier );///读取文件
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
	////图像预处理
	//typedef itk::CastImageFilter<
	//ImageType, FloatImageType >  CastFilterType;
	//通过调用New( )操作来创建对象滤波器并将结果指向itk::SmartPointers

	CastFilterType::Pointer       allT1_castFilter       = CastFilterType::New();
	allT1_castFilter->SetInput(       allT1reader->GetOutput() );///将读取到的数据T1reader通过castfiltertype得到滤波后的图像
	allT1_castFilter->Update();


	
	allT2_castFilter->SetInput(       allT2reader->GetOutput() );///将read的图作为输入
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
	DWI_filter->SetOutputMaximum( 255 );///归到0-255
	DWI_filter->SetInput( allT2_castFilter->GetOutput());
	
	DWI_filter->Update();


	//ITK TO VTK
	typedef itk::ImageToVTKImageFilter< FloatImageType> itkTovtkFilterType;
	itkTovtkFilterType::Pointer itkTovtkImageFilter1 = itkTovtkFilterType::New();
	itkTovtkImageFilter1->SetInput(allT1_castFilter->GetOutput());
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

	//vtkSmartPointer<vtkRenderer> renderer =
	  //  vtkSmartPointer<vtkRenderer>::New();//负责管理场景的渲染过程将场景对象组合起来
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
	renderWindowludi->SetWindowName("颅底配准");
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
	    vtkSmartPointer<vtkRenderWindowInteractor>::New();//提供了一种独立于品台的交互举止相应鼠标案件始终等消息，监听到消息会自动调用回调函数。通过调用invokeevent函数将平台实践翻译成vtk事件
	renderWindowInteractor->SetPicker(pointPicker);
	renderWindowInteractor->SetRenderWindow(renderWindow);

	vtkSmartPointer<PointPickerInteractorStyle> style =
	    vtkSmartPointer<PointPickerInteractorStyle>::New();
	renderWindowInteractor->SetInteractorStyle( style );
	renderWindowInteractor->Start();*/
	
	
}

void saveresulttodicom()//配准后将结果保存成dicom
{
	////////////////////fix
  typedef itk::GDCMImageIO                        ImageIOType;
  typedef itk::GDCMSeriesFileNames                NamesGeneratorType;

  ImageIOType::Pointer gdcmIOf = ImageIOType::New();
  NamesGeneratorType::Pointer namesGeneratorf = NamesGeneratorType::New();
   //namesGeneratorf->SetInputDirectory("C:\\Users\\Administrator\\Desktop\\niuxiahe\\maruohan-data\\haoliangliang(left)\\haoliangliang20151211L\\ZC255IGE\\C0LUTGR3"/*"C:\\Users\\Administrator\\Desktop\\devide\\haolinlin500\\testsquare12"*/ );
  namesGeneratorf->SetInputDirectory(allfixpth);//注意这里是反的
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

	ResampleFilterType::Pointer resampler = ResampleFilterType::New();//变换结果将待配准图像映射到参照图像空间*/

	resamplerdicom->SetTransform( finalTransform );
	resamplerdicom->SetInput( reader->GetOutput() );/////////////全局图像

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
	ImageType > DifferenceFilterType;//减法
	DifferenceFilterType::Pointer difference = DifferenceFilterType::New();

	typedef itk::RescaleIntensityImageFilter<
	ImageType,
	ImageType >   RescalerType;///重新调节亮度以使它们更加明显使负值可视化成为可能
	RescalerType::Pointer intensityRescaler = RescalerType::New();

	intensityRescaler->SetInput( difference->GetOutput() );
	intensityRescaler->SetOutputMinimum(   0 );
	intensityRescaler->SetOutputMaximum( 255 );

	difference->SetInput1( readerf->GetOutput() );//固定图像
	difference->SetInput2( resamplerdicom->GetOutput() );//重新采样后的浮动图像

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
	ImageType >   RescalerType;///重新调节亮度以使它们更加明显使负值可视化成为可能
	RescalerType::Pointer intensityRescalerori = RescalerType::New();

	intensityRescalerori->SetInput( differenceori->GetOutput() );
	intensityRescalerori->SetOutputMinimum(   0 );
	intensityRescalerori->SetOutputMaximum( 255 );

	differenceori->SetInput1( readerf->GetOutput() );//固定图像
	differenceori->SetInput2( reader->GetOutput() );//重新采样后的浮动图像

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
	move_castFilter->SetInput(       reader->GetOutput() );///将读取到的数据T1reader通过castfiltertype得到滤波后的图像
	move_castFilter->Update();

	CastFilterType::Pointer       fix_castFilter       = CastFilterType::New();
	fix_castFilter->SetInput(       readerf->GetOutput() );///将read的图作为输入
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
	    vtkSmartPointer<vtkRenderWindowInteractor>::New();//提供了一种独立于品台的交互举止相应鼠标案件始终等消息，监听到消息会自动调用回调函数。通过调用invokeevent函数将平台实践翻译成vtk事件
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
	
	///////////变换后的图像
	ResampleFilterType::Pointer resampler = ResampleFilterType::New();
	resampler->SetTransform( finalTransform );
	resampler->SetInput( T1_castFilter->GetOutput() );/////////////全局图像

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
	renderWindowresultcon->SetWindowName("局部配准效果");
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
	resamplermha->SetInput( T1_filter->GetOutput() );/////////////全局图像

	FixedImageType::Pointer fixedImage = DWI_filter->GetOutput();

	resamplermha->SetSize(    fixedImage->GetLargestPossibleRegion().GetSize() );
	resamplermha->SetOutputOrigin(  fixedImage->GetOrigin() );
	resamplermha->SetOutputSpacing( fixedImage->GetSpacing() );
	resamplermha->SetOutputDirection( fixedImage->GetDirection() );
	resamplermha->SetDefaultPixelValue( 100 );

	typedef  unsigned char                                          OutputPixelType;
	typedef itk::Image< OutputPixelType, Dimension >                OutputImageType;
	typedef itk::CastImageFilter< FixedImageType, OutputImageType > CastFilterTypeN;//转化重采样的图像的像素类型到最终writer的像素类型
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