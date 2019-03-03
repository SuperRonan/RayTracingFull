#include <Geometry/Texture.h>
#include <Math/Vectorf.h>
#include <Geometry/Ray.h>
#include <Geometry/Triangle.h>
#include <Geometry/CastedRay.h>
#include <stdlib.h>
#include <iostream>
#include <Geometry/RGBColor.h>
#include <Geometry/Material.h>
#include <Geometry/PointLight.h>
#include <Geometry/Camera.h>
#include <Geometry/Cube.h>
#include <Geometry/Disk.h>
#include <Geometry/Cylinder.h>
#include <Geometry/Cone.h>
#include <Visualizer/Visualizer.h>
#include <Geometry/Scene.h>
#include <Geometry/Cornel.h>
#include <Geometry/Loader3ds.h>
#include <Geometry/BoundingBox.h>
#include <omp.h>
#include <Geometry/grid_scene.h>
#include <Geometry/bvh_scene.h>
#include <Geometry/turbo_bvh_scene.h>
#include <chrono>
#include <array>

/// <summary>
/// The directory of the 3D objetcs
/// </summary>
const std::string m_modelDirectory = "..\\..\\Models";

using Geometry::RGBColor;

////////////////////////////////////////////////////////////////////////////////////////////////////
/// \fn	void createGround(Geometry::Scene & scene)
///
/// \brief	Adds a ground to the scene.
///
/// \author	F. Lamarche, Université de Rennes 1
/// \date	03/12/2013
///
/// \param [in,out]	scene	The scene.
////////////////////////////////////////////////////////////////////////////////////////////////////
void createGround(Geometry::Scene & scene)
{
	Geometry::BoundingBox sb = scene.getBoundingBox();
	// Non emissive 
	Geometry::Material * material = new Geometry::Material(RGBColor(), RGBColor(0.3, 0.3, 0.3), RGBColor(0.7, 0.7, 0.7) , 1000.0f); // Non existing material...
	//Geometry::Material * material = new Geometry::Material(RGBColor(), RGBColor(0., 0., 0.), RGBColor(1., 1., 1.), 10000.0f);	//perfect mirror
	//Geometry::Material * material = new Geometry::Material(RGBColor(), RGBColor(1.0f, 1.0f, 1.0f), RGBColor(0.f, 0.f, 0.f), 100.0f); //not a mirror

	//Geometry::Material * material = new Geometry::Material(RGBColor(), RGBColor(1.0,.4,.4)*0.5, RGBColor(), 1000.0f, RGBColor(0.5, 0.5, 0.5)*0.5); // Non existing material...

	Geometry::Square square(material);
	Math::Vector3f scaleV = (sb.max() - sb.min());
	double scale = ::std::max(scaleV[0], scaleV[1])*2.0;
	square.scaleX(scale);
	square.scaleY(scale);
	Math::Vector3f center = (sb.min() + sb.max()) / 2.0;
	center[2] = sb.min()[2];
	square.translate(center);
	scene.add(square);
	::std::cout << "Bounding box: " << sb.min() << "/ " << sb.max() << ", scale: " << scale << ::std::endl;
	::std::cout << "center: " << (sb.min() + sb.max()) / 2.0 << ::std::endl;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// \fn	void createGround(Geometry::Scene & scene)
///
/// \brief	Adds a sirface area light to the scene
///
/// \author	F. Lamarche, Université de Rennes 1
/// \date	03/12/2013
///
/// \param [in,out]	scene	The scene.
////////////////////////////////////////////////////////////////////////////////////////////////////
void createSurfaceLigth(Geometry::Scene & scene, double value)
{
	Geometry::BoundingBox sb = scene.getBoundingBox();
	Geometry::Material * material = new Geometry::Material(RGBColor(), RGBColor(), RGBColor(), 100.0f, RGBColor(value, value, value));
	Geometry::Square square(material);
	Math::Vector3f scaleV = (sb.max() - sb.min());
	//double scale = ::std::max(scaleV[0], scaleV[1])*0.1;
	double factor = 0.5;
	square.scaleX(scaleV[0] * factor);
	square.scaleY(scaleV[1] * factor);
	Math::Vector3f center = (sb.min() + sb.max()) / 2.0;
	center[2] = sb.max()[2] * 3;
	square.translate(center);
	scene.add(square);
}






void init_test(Geometry::Scene & scene, Visualizer::Visualizer const& visu)
{
	if(true)
	{
		Geometry::Material * material = new Geometry::Material(RGBColor(), RGBColor(1, 1, 1));

		Geometry::OneTriangle t(material);
		t.scale(10);
		t.translate(Math::makeVector(10, 0, 0));
		scene.add(t);

		Geometry::OneTriangle t2(material);
		t2.scale(5);
		t2.rotate(Math::Quaternion<double>(Math::makeVector(0, 1, 0), Math::piDiv2));
		t2.translate(Math::makeVector(0.0, 0.0, 2.5));

		scene.add(t2);

		Geometry::Square square(material);
		square.scale(10);

		scene.add(square);



		// 2.2 Adds point lights in the scene 
		{
			Geometry::PointLight pointLight(Math::makeVector(5.0f, 0.f, 2.0f), RGBColor(1.f, 1.f, 1.f));
			scene.add(pointLight);
		}
		{
			Geometry::Camera camera(Math::makeVector(5.f, 0.1f, 2.0f), Math::makeVector(0.0f, 0.0f, 0.0f), 0.3f, ((double)visu.width()) / ((double)visu.height()), 1.0f);
			scene.setCamera(camera);
		}
	}
	else
	{

	}
}








////////////////////////////////////////////////////////////////////////////////////////////////////
/// \fn	void initDiffuse(Geometry::Scene & scene)
///
/// \brief	Adds a Cornell Box with diffuse material on each walls to the scene. This Cornel box
/// 		contains two cubes.
///
/// \author	F. Lamarche, Université de Rennes 1
/// \date	03/12/2013
///
/// \param [in,out]	scene	The scene.
////////////////////////////////////////////////////////////////////////////////////////////////////
void initDiffuse(Geometry::Scene & scene, Visualizer::Visualizer const& visu)
{
	Geometry::Material * material = new Geometry::Material(RGBColor(), RGBColor(0, 0, 0.0), RGBColor(0.95f, 0.95f, 0.95f), 1, RGBColor());
	Geometry::Material * material2 = new Geometry::Material(RGBColor(), RGBColor(1.0, 1.0, 1.0), RGBColor(0, 0, 0), 1000, RGBColor());
	Geometry::Material * cubeMat = new Geometry::Material(RGBColor(), RGBColor(1, 0.0, 0.0), RGBColor(0.0, 0.0, 0.0), 20.0f, RGBColor());
	//Geometry::Material * cubeMat = new Geometry::Material(RGBColor(), RGBColor(1.0f,0.0,0.0), RGBColor(0.0,0.0,0.0), 20.0f, RGBColor(10.0,0,0)) ;
	Geometry::Cornel geo(material2, material2, material, material, material, material);

	geo.scaleX(10);
	geo.scaleY(10);
	geo.scaleZ(10);
	scene.add(geo);

	Geometry::Cube tmp(cubeMat);
	tmp.translate(Math::makeVector(1.5, -1.5, 0.0));
	scene.add(tmp);

	Geometry::Cube tmp2(cubeMat);
	tmp2.translate(Math::makeVector(2, 1, -4));
	scene.add(tmp2);

	

	// 2.2 Adds point lights in the scene 
	{
		Geometry::PointLight pointLight(Math::makeVector(0.0f, 0.f, 2.0f), RGBColor(0.5f, 0.5f, 0.5f));
		scene.add(pointLight);
	}
	{
		Geometry::PointLight pointLight2(Math::makeVector(4.f, 0.f, 0.f), RGBColor(0.5f, 0.5f, 0.5f));
		scene.add(pointLight2);
	}
	{
		Geometry::Camera camera(Math::makeVector(-4.0f, 0.0f, 0.0f), Math::makeVector(0.0f, 0.0f, 0.0f), 0.3f, ((double)visu.width())/((double)visu.height()), 1.0f);
		scene.setCamera(camera);
	}
}




void initDiffuse_surface_light(Geometry::Scene & scene, Visualizer::Visualizer const& visu)
{
	Geometry::Material * material = new Geometry::Material(RGBColor(), RGBColor(0, 0, 0.0), RGBColor(0.95f, 0.95f, 0.95f), 1, RGBColor());
	Geometry::Material * material2 = new Geometry::Material(RGBColor(), RGBColor(1.0, 1.0, 1.0), RGBColor(0, 0, 0), 1000, RGBColor());
	Geometry::Material * cubeMat = new Geometry::Material(RGBColor(), RGBColor(10, 0.0, 0.0), RGBColor(0.0, 0.0, 0.0), 20.0f, RGBColor());
	//Geometry::Material * cubeMat = new Geometry::Material(RGBColor(), RGBColor(1.0f,0.0,0.0), RGBColor(0.0,0.0,0.0), 20.0f, RGBColor(10.0,0,0)) ;
	Geometry::Cornel geo(material2, material2, material2, material2, material2, material2);

	geo.scaleX(10);
	geo.scaleY(10);
	geo.scaleZ(10);
	scene.add(geo);

	Geometry::Cube tmp(cubeMat);
	tmp.scale(2);
	tmp.translate(Math::makeVector(3.5, -3.5, 0.0));
	scene.add(tmp);

	Geometry::Cube tmp2(cubeMat);
	tmp2.translate(Math::makeVector(2, 1, -4));
	scene.add(tmp2);



	// 2.2 Adds point lights in the scene 
	{
		Geometry::PointLight pointLight(Math::makeVector(0.0f, 0.f, 2.0f), RGBColor(0.5f, 0.5f, 0.5f));
		scene.add(pointLight);
	}
	{
		Geometry::PointLight pointLight2(Math::makeVector(4.f, 0.f, 0.f), RGBColor(0.5f, 0.5f, 0.5f));
		scene.add(pointLight2);
	}
	//Add the surface light
	{
		Geometry::Material * surface = new Geometry::Material(0, 0, 0, 0, { 1, 1, 1 });
		Geometry::Square slight(surface);
		slight.scaleY(2);
		slight.scaleX(2);
		slight.translate(Math::makeVector(0.0, 0.0, 4.5));
		scene.add(slight);
	}

	/*
	{
		Geometry::Material * surface = new Geometry::Material(0, 0, 0, 0, { 1, 1, 1 });
		Geometry::OneTriangle slight(surface);
		slight.scaleY(4.0);
		slight.scaleX(4.0);
		slight.translate(Math::makeVector(0.0, 2.0, -4.9));
		scene.add(slight);
	}
	*/


	//Add the Camera
	{
		Geometry::Camera camera(Math::makeVector(-4.0f, 0.0f, 0.0f), Math::makeVector(0.0f, 0.0f, 0.0f), 0.3f, ((double)visu.width()) / ((double)visu.height()), 1.0f);
		scene.setCamera(camera);
	}

	
}




void initDiffuse_dif(Geometry::Scene & scene, Visualizer::Visualizer const& visu)
{
	Geometry::Material * material = new Geometry::Material(RGBColor(), RGBColor(0, 0, 0.0), RGBColor(0.95f, 0.95f, 0.95f), 1, RGBColor());
	Geometry::Material * material2 = new Geometry::Material(RGBColor(), RGBColor(1.0, 1.0, 1.0), RGBColor(0, 0, 0), 1000, RGBColor());
	Geometry::Material * cubeMat = new Geometry::Material(RGBColor(), RGBColor(1, 0.0, 0.0), RGBColor(0.0, 0.0, 0.0), 20.0f, RGBColor());
	//Geometry::Material * cubeMat = new Geometry::Material(RGBColor(), RGBColor(1.0f,0.0,0.0), RGBColor(0.0,0.0,0.0), 20.0f, RGBColor(10.0,0,0)) ;
	Geometry::Cornel geo(material2, material2, material, material, material, material);

	geo.scaleX(10);
	geo.scaleY(10);
	geo.scaleZ(10);
	scene.add(geo);

	Geometry::Cube tmp(cubeMat);
	tmp.translate(Math::makeVector(1.5, -1.5, 0.0));
	scene.add(tmp);

	Geometry::Cube tmp2(cubeMat);
	tmp2.translate(Math::makeVector(2, 1, -4));
	scene.add(tmp2);

	Geometry::Material * transparent_mat = new Geometry::Material(RGBColor(), RGBColor(), RGBColor(0.95, 0.95, 0.95), 1000);
	transparent_mat->m_transparent = true;

	Geometry::Cylinder transparent_obj(20, 0.5, 0.5, transparent_mat);
	transparent_obj.scale(1.5);
	transparent_obj.m_medium = Geometry::medium(2.0);
	transparent_mat->m_medium = transparent_obj.get_medium();
	scene.add(transparent_obj);


	// 2.2 Adds point lights in the scene 
	{
		Geometry::PointLight pointLight(Math::makeVector(-4.0f, 0.f, 4.0f), RGBColor(0.5f, 0.5f, 0.5f));
		scene.add(pointLight);
	}
	{
		Geometry::PointLight pointLight2(Math::makeVector(4.f, 0.f, 0.f), RGBColor(0.5f, 0.5f, 0.5f));
		scene.add(pointLight2);
	}
	{
		Geometry::Camera camera(Math::makeVector(-4.0f, 0.0f, 0.0f), Math::makeVector(0.0f, 0.0f, 0.0f), 0.3f, ((double)visu.width()) / ((double)visu.height()), 1.0f);
		scene.setCamera(camera);
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// \fn	void initSpecular(Geometry::Scene & scene)
///
/// \brief	Adds a Cornel box in the provided scene. Walls are specular and the box contains two 
/// 		cubes.
///
/// \author	F. Lamarche, Université de Rennes 1
/// \date	03/12/2013
////////////////////////////////////////////////////////////////////////////////////////////////////
void initSpecular(Geometry::Scene & scene, Visualizer::Visualizer const& visu)
{
	Geometry::Material * material = new Geometry::Material(RGBColor(), RGBColor(0, 0, 0.0), RGBColor(0.7f, 0.7f, 0.7f), 100, RGBColor());
	Geometry::Material * material2 = new Geometry::Material(RGBColor(), RGBColor(0, 0, 1.0f), RGBColor(0, 0, 0), 1000, RGBColor());
	//Geometry::Material * cubeMat = new Geometry::Material(RGBColor(), RGBColor(1.0f,0.0,0.0), RGBColor(0.0,0.0,0.0), 20.0f, RGBColor(10.0,0,0)) ;
	Geometry::Material * cubeMat = new Geometry::Material(RGBColor(), RGBColor(1.0f, 0.0, 0.0), RGBColor(0.0, 0.0, 0.0), 20.0f, RGBColor());
	Geometry::Cornel geo(material, material, material, material, material, material); //new Geometry::Cube(material2) ;////new Cone(4, material) ; //new Geometry::Cylinder(5, 1, 1, material) ;////////new Geometry::Cube(material) ;////; //new Geometry::Cube(material) ; //new Geometry::Cylinder(100, 2, 1, material) ; //

	

	geo.scaleX(10);
	geo.scaleY(10);
	geo.scaleZ(10);
	scene.add(geo);

	Geometry::Cube tmp(cubeMat);
	tmp.translate(Math::makeVector(1.5, -1.5, 0.0));
	scene.add(tmp);

	Geometry::Cube tmp2(cubeMat);
	tmp2.translate(Math::makeVector(2, 1, -4));
	scene.add(tmp2);

	// 2.2 Adds point lights in the scene 
	{
		Geometry::PointLight pointLight(Math::makeVector(0.0f, 0.f, 2.0f), RGBColor(0.5f, 0.5f, 0.5f) * 5);
		scene.add(pointLight);
	}
	{
		Geometry::PointLight pointLight2(Math::makeVector(4.f, 0.f, 0.f), RGBColor(0.5f, 0.5f, 0.5f) * 5);
		scene.add(pointLight2);
	}
	// Sets the camera
	{
		Geometry::Camera camera(Math::makeVector(-4.0f, 0.0f, 0.0f), Math::makeVector(0.0f, 0.0f, 0.0f), 0.3f, ((double)visu.width()) / ((double)visu.height()), 1.0f);
		scene.setCamera(camera);
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// \fn	void initDiffuseSpecular(Geometry::Scene & scene)
///
/// \brief	Adds a Cornel box in the provided scene. The cornel box as diffuse and specular walls and 
/// 		contains two boxes.
///
/// \author	F. Lamarche, Université de Rennes 1
/// \date	03/12/2013
///
/// \param [in,out]	scene	The scene.
////////////////////////////////////////////////////////////////////////////////////////////////////
void initDiffuseSpecular(Geometry::Scene & scene, Visualizer::Visualizer const& visu)
{
	Geometry::Material * material = new Geometry::Material(RGBColor(), RGBColor(0, 0, 0.0), RGBColor(0.7f, 0.71f, 0.7f), 1000, RGBColor());
	Geometry::Material * material2 = new Geometry::Material(RGBColor(), RGBColor(1, 1, 1.0f), RGBColor(0, 0, 0), 1000, RGBColor());
	Geometry::Material * cubeMat = new Geometry::Material(RGBColor(), RGBColor(0., 0., 0.), RGBColor(0.7, 0.4, 0.), 20.0f, RGBColor());//gold
	Geometry::Material * cubeMat3 = new Geometry::Material(RGBColor(), RGBColor(0., 0., 0.), RGBColor(0.5, 0.5, 0.5), 20.0f, RGBColor());//silver
	Geometry::Material * cubeMat2 = new Geometry::Material(RGBColor(), RGBColor(1.0f, 0.0, 0.0), RGBColor(0.0, 0.0, 0.0), 20.0f, RGBColor());//red
																																			 //Geometry::Material * cubeMat = new Geometry::Material(RGBColor(), RGBColor(0.0f,0.0,0.0), RGBColor(0.0,0.0,0.0), 20.0f, RGBColor(10.0,0,0)) ;
																																			 //Geometry::Material * cubeMat2 = new Geometry::Material(RGBColor(), RGBColor(0.0f,0.0,0.0), RGBColor(0.0,0.0,0.0), 20.0f, RGBColor(0.0,10,0)) ;
	Geometry::Cornel geo(material2, material2, material, material, material, material); //new Geometry::Cube(material2) ;////new Cone(4, material) ; //new Geometry::Cylinder(5, 1, 1, material) ;////////new Geometry::Cube(material) ;////; //new Geometry::Cube(material) ; //new Geometry::Cylinder(100, 2, 1, material) ; //

	Geometry::Cylinder cyl(100, 1., -1., cubeMat3);
	cyl.translate(Math::makeVector(0., -2., -2.));
	
	scene.add(cyl);

	geo.scaleX(10);
	geo.scaleY(10);
	geo.scaleZ(10);
	scene.add(geo);


	Geometry::Cube tmp(cubeMat2);
	tmp.translate(Math::makeVector(1.5, -1.5, 0.0));
	scene.add(tmp);

	Geometry::Cube tmp2(cubeMat);
	tmp2.scaleY(5);
	tmp2.scaleX(4);
	tmp2.translate(Math::makeVector(2, 1, -4));

	scene.add(tmp2);

	// 2.2 Adds point lights in the scene 
	{
		Geometry::PointLight pointLight(Math::makeVector(0.0f, 0.f, 2.0f), RGBColor(0.5f, 0.5f, 0.5f)*5);
		scene.add(pointLight);
	}
	{
		Geometry::PointLight pointLight2(Math::makeVector(4.f, 0.f, 0.f), RGBColor(0.5f, 0.5f, 0.5f)*5);
		scene.add(pointLight2);
	}
	{
		Geometry::Camera camera(Math::makeVector(-4.f, .0f, .0f), Math::makeVector(0.0f, 0.0f, 0.0f), 0.4f, ((double)visu.width()) / ((double)visu.height()), 1.0f);
		scene.setCamera(camera);
	}


	Geometry::Material * black = new Geometry::Material(RGBColor(), RGBColor(0.1, 0.1, 0.1), RGBColor(0., 0., 0.), 10, RGBColor());
	{
		Geometry::Cube pillar(black);

		pillar.scaleX(.1);
		pillar.scaleY(.1);
		pillar.scaleZ(10);

		pillar.translate(Math::makeVector(5., 5., 0.));
		scene.add(pillar);

		pillar.translate(Math::makeVector(0., -10., 0.));
		scene.add(pillar);

		pillar.translate(Math::makeVector(-10., 0., 0.));
		scene.add(pillar);

		pillar.translate(Math::makeVector(0., 10., 0.));
		scene.add(pillar);
	}
	{
		Geometry::Cube pillar(black);

		pillar.scaleX(10);
		pillar.scaleY(.1);
		pillar.scaleZ(.1);

		pillar.translate(Math::makeVector(0., 5., 5.));
		scene.add(pillar);

		pillar.translate(Math::makeVector(0., 0., -10.));
		scene.add(pillar);


		pillar.translate(Math::makeVector(0., -10., 0.));
		scene.add(pillar);

		pillar.translate(Math::makeVector(0., 0., 10.));
		scene.add(pillar);
	}
	{
		Geometry::Cube pillar(black);

		pillar.scaleX(0.1);
		pillar.scaleY(10);
		pillar.scaleZ(.1);

		pillar.translate(Math::makeVector(5., 0., 5.));
		scene.add(pillar);

		pillar.translate(Math::makeVector(0., 0., -10.));
		scene.add(pillar);


		pillar.translate(Math::makeVector(-10., 0., 0.));
		scene.add(pillar);

		pillar.translate(Math::makeVector(0., 0., 10.));
		scene.add(pillar);
	}
}

/// <summary>
/// Intializes a scene containing a garage.
/// </summary>
/// <param name="scene"></param>
void initGarage(Geometry::Scene & scene, Visualizer::Visualizer const& visu)
{
	Geometry::Loader3ds loader(m_modelDirectory + "\\garage.3ds", m_modelDirectory + "");

	for (size_t cpt = 0; cpt < loader.getMeshes().size(); ++cpt)
	{
		scene.add(*loader.getMeshes()[cpt]);
	}

	// 2.2 Adds point lights in the scene 
	Geometry::BoundingBox sb = scene.getBoundingBox();
	Math::Vector3f position = sb.max();
	Geometry::PointLight light1(position, RGBColor(1.0, 1.0, 1.0) * 500);
	scene.add(light1);
	position[1] = -position[1];
	Geometry::PointLight light2(position, RGBColor(1.0, 1.0, 1.0) * 500);
	scene.add(light2);
	{
		Geometry::Camera camera(Math::makeVector(750.0f, -1500.f, 1000.f)*0.85f, Math::makeVector(200.0f, 0.0f, 0.0f), 0.3f, ((double)visu.width()) / ((double)visu.height()), 1.0f);
		scene.setCamera(camera);
	}
	createGround(scene);
}

/// <summary>
/// Initializes a scene containing a guitar.
/// </summary>
/// <param name="scene"></param>
void initGuitar(Geometry::Scene & scene, Visualizer::Visualizer const& visu)
{
	Geometry::Loader3ds loader(m_modelDirectory + "\\guitar2.3ds", m_modelDirectory + "");

	for (size_t cpt = 0; cpt < loader.getMeshes().size(); ++cpt)
	{
		scene.add(*loader.getMeshes()[cpt]);
	}
	// 2.2 Adds point lights in the scene 
	Geometry::BoundingBox sb = scene.getBoundingBox();
	Math::Vector3f position = sb.max();
	Geometry::PointLight light1(position + Math::makeVector(0.f, 0.f, 70.f/*100.f*/), RGBColor(1.0, 1.0, 1.0)*200 );
	scene.add(light1);
	position[1] = -position[1];
	Geometry::PointLight light2(position + Math::makeVector(0.f, 0.f, 200.f), RGBColor(1.0, 1.0, 1.0)*200 );
	scene.add(light2);
	{
		Geometry::Camera camera(Math::makeVector(-500., -1000., 1000.)*1.05, Math::makeVector(500.f, 0.0f, 0.0f), 0.6f, ((double)visu.width()) / ((double)visu.height()), 1.0f);
		camera.translateLocal(Math::makeVector(100.0f, -100.f, -200.f));
		scene.setCamera(camera);
	}
	createGround(scene);
}

/// <summary>
/// Initializes a scene containing a dog
/// </summary>
/// <param name="scene"></param>
void initDog(Geometry::Scene & scene, Visualizer::Visualizer const& visu)
{
	Geometry::Loader3ds loader(m_modelDirectory + "\\Dog\\dog.3ds", m_modelDirectory + "\\dog");

	for (size_t cpt = 0; cpt < loader.getMeshes().size(); ++cpt)
	{
		scene.add(*loader.getMeshes()[cpt]);
	}

	// 2.2 Adds point lights in the scene 
	Geometry::BoundingBox sb = scene.getBoundingBox();
	Math::Vector3f position = sb.max();
	Geometry::PointLight light1(position, RGBColor(1.0, 1.0, 1.0) *.5);
	scene.add(light1);
	position[1] = -position[1];
	Geometry::PointLight light2(position, RGBColor(1.0, 1.0, 1.0) *.5);
	scene.add(light2);
	{
		Geometry::Camera camera(Math::makeVector(10.f, 10.f, 6.f)*0.5, Math::makeVector(0.f, 0.0f, 2.5f), .7f, ((double)visu.width()) / ((double)visu.height()), 1.0f);
		camera.translateLocal(Math::makeVector(-0.4f, 0.f, 0.9f));
		scene.setCamera(camera);
	}
	createGround(scene);
}

/// <summary>
/// Initializes a scene containing a temple.
/// </summary>
/// <param name="scene"></param>
void initTemple(Geometry::Scene & scene, Visualizer::Visualizer const& visu)
{
	Geometry::BoundingBox box(Math::makeVector(0.0f, 0.0f, 0.0f), Math::makeVector(0.0f, 0.0f, 0.0f));
	Geometry::Loader3ds loader(m_modelDirectory + "\\Temple\\Temple of St Seraphim of Sarov N270116_2.3ds", m_modelDirectory + "\\Temple");

	for (size_t cpt = 0; cpt < loader.getMeshes().size(); ++cpt)
	{
		//loader.getMeshes()[cpt]->translate(Math::makeVector(20.f, 0.f, 40.0f));
		scene.add(*loader.getMeshes()[cpt]);
		box.update(*loader.getMeshes()[cpt]);
	}


	// 2.2 Adds point lights in the scene 
	Geometry::BoundingBox sb = scene.getBoundingBox();
	Math::Vector3f position = sb.max();
	Geometry::PointLight light1(position, RGBColor(1.0, 1.0, 1.0)*20.0);
	//Geometry::PointLight light1(position, RGBColor(1.0, 1.0, 1.0));
	scene.add(light1);
	position[1] = -position[1];
	Geometry::PointLight light2(position, RGBColor(1.0, 1.0, 1.0)*10);
	//Geometry::PointLight light2(position, RGBColor(1.0, 1.0, 1.0) );
	scene.add(light2);
	{
		Geometry::Camera camera(Math::makeVector(20.0f, -100.0f, 15.0f), Math::makeVector(-20.f, 0.f, -40.f), 0.3f, ((double)visu.width()) / ((double)visu.height()), 1.0f);
		camera.translateLocal(Math::makeVector(40.f, 0.f, 0.f));
		scene.setCamera(camera);
	}
	createGround(scene);
	createSurfaceLigth(scene, 50);
	//createSurfaceLigth(scene, 500);
}

/// <summary>
/// Initializes a scene containing a robot.
/// </summary>
/// <param name="scene"></param>
void initRobot(Geometry::Scene & scene, Visualizer::Visualizer const& visu)
{
	Geometry::BoundingBox box(Math::makeVector(0.0f, 0.0f, 0.0f), Math::makeVector(0.0f, 0.0f, 0.0f));
	Geometry::Loader3ds loader(m_modelDirectory + "\\Robot.3ds", m_modelDirectory + "");

	for (size_t cpt = 0; cpt < loader.getMeshes().size(); ++cpt)
	{
		//loader.getMeshes()[cpt]->translate(Math::makeVector(20.f, 0.f, 40.0f));
		scene.add(*loader.getMeshes()[cpt]);
		box.update(*loader.getMeshes()[cpt]);
	}

	// 2.2 Adds point lights in the scene 
	Geometry::BoundingBox sb = scene.getBoundingBox();
	Math::Vector3f position = sb.max();
	Geometry::PointLight light1(position, RGBColor(1.0, 1.0, 1.0) * 10);
	scene.add(light1);
	position[1] = -position[1];
	Geometry::PointLight light2(position, RGBColor(1.0, 1.0, 1.0) *20);
	scene.add(light2);
	{
		Geometry::Camera camera(Math::makeVector(100.0f, -50.0f, 0.0f), Math::makeVector(0.f, 0.f, 5.f), 0.3f, ((double)visu.width()) / ((double)visu.height()), 1.0f);
		camera.translateLocal(Math::makeVector(0.f, 40.f, 50.f));
		scene.setCamera(camera);
	}
	createGround(scene);
}

/// <summary>
/// Initializes a scene containing a grave stone
/// </summary>
/// <param name="scene"></param>
void initGraveStone(Geometry::Scene & scene, Visualizer::Visualizer const& visu)
{
	Geometry::BoundingBox box(Math::makeVector(0.0f, 0.0f, 0.0f), Math::makeVector(0.0f, 0.0f, 0.0f));
	Geometry::Loader3ds loader(m_modelDirectory + "\\gravestone\\GraveStone.3ds", m_modelDirectory + "\\gravestone");

	::std::vector<Geometry::Material*> materials = loader.getMaterials();
	for (auto it = materials.begin(), end = materials.end(); it != end; ++it)
	{
		//(*it)->setAmbient(RGBColor());
		(*it)->setSpecular(RGBColor());
		//(*it)->setSpecular((*it)->getSpecular()*0.05);
	}

	for (size_t cpt = 0; cpt < loader.getMeshes().size(); ++cpt)
	{
		//loader.getMeshes()[cpt]->translate(Math::makeVector(20.f, 0.f, 40.0f));
		scene.add(*loader.getMeshes()[cpt]);
		box.update(*loader.getMeshes()[cpt]);
	}

	// 2.2 Adds point lights in the scene 
	Geometry::BoundingBox sb = scene.getBoundingBox();
	Math::Vector3f position = sb.max();
	Geometry::PointLight light1(position, RGBColor(1.0, 1.0, 1.0)*100 );
	scene.add(light1);
	position[1] = -position[1];
	Geometry::PointLight light2(position, RGBColor(1.0, 1.0, 1.0)*100 );
	scene.add(light2);
	{
		Geometry::Camera camera(Math::makeVector(0.f, -300.0f, 200.0f), Math::makeVector(0.f, 0.f, 125.f), 0.3f, ((double)visu.width()) / ((double)visu.height()), 1.0f);
		camera.translateLocal(Math::makeVector(0.f, 80.f, 120.f));
		scene.setCamera(camera);
	}
	createGround(scene);
	//createSurfaceLigth(scene, 400);
}

/// <summary>
/// Initializes a scene containing a boat
/// </summary>
/// <param name="scene"></param>
void initBoat(Geometry::Scene & scene, Visualizer::Visualizer const& visu)
{
	Geometry::BoundingBox box(Math::makeVector(0.0f, 0.0f, 0.0f), Math::makeVector(0.0f, 0.0f, 0.0f));
	Geometry::Loader3ds loader(m_modelDirectory + "\\Boat\\boat.3ds", m_modelDirectory + "\\boat");
	//// We remove the specular components of the materials...
	::std::vector<Geometry::Material*> materials = loader.getMaterials();
	for (auto it = materials.begin(), end = materials.end(); it != end; ++it)
	{
		//(*it)->setAmbient(RGBColor(0.01, 0.01, 0.01));
		(*it)->setSpecular(RGBColor());
	}

	for (size_t cpt = 0; cpt < loader.getMeshes().size(); ++cpt)
	{
		//loader.getMeshes()[cpt]->translate(Math::makeVector(20.f, 0.f, 40.0f));
		scene.add(*loader.getMeshes()[cpt]);
		box.update(*loader.getMeshes()[cpt]);
	}

	// 2.2 Adds point lights in the scene 
	Geometry::BoundingBox sb = scene.getBoundingBox();
	Math::Vector3f position = sb.max();
	Geometry::PointLight light1(position, RGBColor(1.0, 1.0, 1.0)*6000);
	scene.add(light1);
	position[1] = -position[1];
	Geometry::PointLight light2(position, RGBColor(1.0, 1.0, 1.0) * 20000);
	scene.add(light2);
	{
		Geometry::Camera camera(Math::makeVector(5000.f, 5000.f, 200.0f), Math::makeVector(0.f, 0.f, 60.f), 0.3f, ((double)visu.width()) / ((double)visu.height()), 1.0f);
		camera.translateLocal(Math::makeVector(2000.f, 3000.f, 700.f));
		scene.setCamera(camera);
	}
	createGround(scene);
}

/// <summary>
/// Initializes a scene containing a tibet house
/// </summary>
/// <param name="scene"></param>
void initTibetHouse(Geometry::Scene & scene, Visualizer::Visualizer const& visu)
{
	Geometry::BoundingBox box(Math::makeVector(0.0f, 0.0f, 0.0f), Math::makeVector(0.0f, 0.0f, 0.0f));
	Geometry::Loader3ds loader(m_modelDirectory + "\\TibetHouse\\TibetHouse.3ds", m_modelDirectory + "\\TibetHouse");
	// We remove the specular components of the materials... and add surface light sources :)
	::std::vector<Geometry::Material*> materials = loader.getMaterials();
	for (auto it = materials.begin(), end = materials.end(); it != end; ++it)
	{
		(*it)->setSpecular(RGBColor());
		(*it)->setAmbient(RGBColor());
		//std::cout << (*it)->getAmbient() << std::endl;
		if ((*it)->getTextureFile() == m_modelDirectory + "\\TibetHouse" + "\\3D69C2DE.png")
		{
			(*it)->setEmissive(RGBColor(1.0, 1.0, 1.0)*1);
		}
	}

	for (size_t cpt = 0; cpt < loader.getMeshes().size(); ++cpt)
	{
		//loader.getMeshes()[cpt]->translate(Math::makeVector(20.f, 0.f, 40.0f));
		scene.add(*loader.getMeshes()[cpt]);
		box.update(*loader.getMeshes()[cpt]);
	}

	// 2.2 Adds point lights in the scene 
	Geometry::BoundingBox sb = scene.getBoundingBox();
	Math::Vector3f position = sb.max();
	Geometry::PointLight light1(position, RGBColor(1.0, 1.0, 1.0) * 20);
	scene.add(light1);
	position[1] = -position[1];
	Geometry::PointLight light2(position, RGBColor(1.0, 1.0, 1.0) *20);
	scene.add(light2);
	Geometry::PointLight light3(Math::makeVector(5.f, 35.f, 5.f), RGBColor(1.0, 1.0, 1.0)*20 ); //*50
	scene.add(light3);
	{
		Geometry::Camera camera(Math::makeVector(20.f, 0.f, 0.0f), Math::makeVector(5.f, 35.f, 0.f), 0.3f, ((double)visu.width()) / ((double)visu.height()), 1.0f);
		camera.translateLocal(Math::makeVector(0.f, 5.f, 0.f));// +Math::makeVector(0.0, 5.0, 0.0));
		scene.setCamera(camera);
	}
	createGround(scene);
}

/// <summary>
/// Initializes a scene containing a tibet house. Camera is placed inside the house.
/// </summary>
/// <param name="scene"></param>
void initTibetHouseInside(Geometry::Scene & scene, Visualizer::Visualizer const& visu)
{
	Geometry::BoundingBox box(Math::makeVector(0.0f, 0.0f, 0.0f), Math::makeVector(0.0f, 0.0f, 0.0f));
	Geometry::Loader3ds loader(m_modelDirectory + "\\TibetHouse\\TibetHouse.3ds", m_modelDirectory + "\\TibetHouse");
	// We remove the specular components of the materials... and add surface light sources :)
	::std::vector<Geometry::Material*> materials = loader.getMaterials();
	for (auto it = materials.begin(), end = materials.end(); it != end; ++it)
	{
		(*it)->setSpecular(RGBColor());
		//(*it)->setAmbient(RGBColor());
		if ((*it)->getTextureFile() == m_modelDirectory + "\\TibetHouse" + "\\3D69C2DE.png")
		{
			(*it)->setEmissive(RGBColor(1.0, .2, 0.0));
		}
	}

	for (size_t cpt = 0; cpt < loader.getMeshes().size(); ++cpt)
	{
		//loader.getMeshes()[cpt]->translate(Math::makeVector(20.f, 0.f, 40.0f));
		scene.add(*loader.getMeshes()[cpt]);
		box.update(*loader.getMeshes()[cpt]);
	}

	// 2.2 Adds point lights in the scene 
	Geometry::BoundingBox sb = scene.getBoundingBox();
	Math::Vector3f position = sb.max();
	Geometry::PointLight light1(position, RGBColor(1.0, 1.0, 1.0)*10);
	scene.add(light1);
	position[1] = -position[1];
	Geometry::PointLight light2(position, RGBColor(1.0, 1.0, 1.0)*10);
	scene.add(light2);
	Geometry::PointLight light3(Math::makeVector(5.f, 35.f, 5.f), RGBColor(1.0, 1.0, 1.0)*10 );
	scene.add(light3);
	{
		Geometry::Camera camera(Math::makeVector(20.f, 0.f, 5.0f), Math::makeVector(5.f, 35.f, 5.f), 0.3f, ((double)visu.width()) / ((double)visu.height()), 1.0f);
		camera.translateLocal(Math::makeVector(10.f - 2, 30.f, 0.f));
		scene.setCamera(camera);
	}
	createGround(scene);
}

/// <summary>
/// Initializes a scene containing a medieval city
/// </summary>
/// <param name="scene"></param>
void initMedievalCity(Geometry::Scene & scene, Visualizer::Visualizer const& visu)
{
	Geometry::BoundingBox box(Math::makeVector(0.0f, 0.0f, 0.0f), Math::makeVector(0.0f, 0.0f, 0.0f));
	Geometry::Loader3ds loader(m_modelDirectory + "\\Medieval\\MedievalCity.3ds", m_modelDirectory + "\\Medieval\\texture");
	// We remove the specular components of the materials...
	::std::vector<Geometry::Material*> materials = loader.getMaterials();
	for (auto it = materials.begin(), end = materials.end(); it != end; ++it)
	{
		(*it)->setSpecular(RGBColor());
		(*it)->setAmbient(RGBColor());

	}

	for (size_t cpt = 0; cpt < loader.getMeshes().size(); ++cpt)
	{
		//loader.getMeshes()[cpt]->translate(Math::makeVector(20.f, 0.f, 40.0f));
		//loader.getMeshes()[cpt]->rotate(Math::Quaternion<double>(Math::makeVector(0.0, 1.0, 0.0), Math::pi));
		scene.add(*loader.getMeshes()[cpt]);
		box.update(*loader.getMeshes()[cpt]);
	}

	// 2.2 Adds point lights in the scene 
	Geometry::BoundingBox sb = scene.getBoundingBox();
	Math::Vector3f position = sb.max();
	Geometry::PointLight light1(position + Math::makeVector(0.0, 0.0, 150.0), RGBColor(1.0, 0.6, 0.3) * 100);
	scene.add(light1);
	position[1] = -position[1];
	Geometry::PointLight light2(position + Math::makeVector(0.0, 0.0, 1000.0), RGBColor(1.0, 1.0, 1.0) * 100);
	//scene.add(light2);
	Geometry::PointLight light3(Math::makeVector(5.f, 35.f, 5.f), RGBColor(1.0, 1.0, 1.0) * 50);
	//scene.add(light3);
	{
		Geometry::Camera camera(Math::makeVector(0.f, 300.f, 1000.0f), Math::makeVector(0.0, 0.0, 0.0), 0.3f, ((double)visu.width()) / ((double)visu.height()), 1.0f);
		camera.translateLocal(Math::makeVector(0.0, 800., -100.0));
		scene.setCamera(camera);
	}
	createGround(scene);
}

/// <summary>
/// Initializes a scene containing a sombrero
/// </summary>
/// <param name="scene"></param>
void initSombrero(Geometry::Scene & scene, Visualizer::Visualizer const& visu)
{
	Geometry::BoundingBox box(Math::makeVector(0.0f, 0.0f, 0.0f), Math::makeVector(0.0f, 0.0f, 0.0f));
	Geometry::Loader3ds loader(m_modelDirectory + "\\sombrero\\sombrero.3ds", m_modelDirectory + "\\sombrero");
	// We remove the specular components of the materials...
	::std::vector<Geometry::Material*> materials = loader.getMaterials();
	for (auto it = materials.begin(), end = materials.end(); it != end; ++it)
	{
		(*it)->setSpecular(RGBColor());
		(*it)->setAmbient(RGBColor());
	}

	for (size_t cpt = 0; cpt < loader.getMeshes().size(); ++cpt)
	{
		scene.add(*loader.getMeshes()[cpt]);
		box.update(*loader.getMeshes()[cpt]);
	}

	// 2.2 Adds point lights in the scene 
	Geometry::BoundingBox sb = scene.getBoundingBox();
	Math::Vector3f position = sb.max();
	Geometry::PointLight light1(position, RGBColor(1.0, 1.0, 1.0)*100);
	scene.add(light1);
	position[1] = -position[1];
	Geometry::PointLight light2(position, RGBColor(1.0, 1.0, 1.0)*100);
	

	Geometry::Material * emisive = new Geometry::Material(0, 0, 0, 0, { 100, 100, 100 });
	Geometry::Square surface(emisive);
	surface.scale(500.0);
	surface.translate(Math::makeVector(0, 0, 400));
	scene.add(surface);
	scene.add(light2);
	{
		Geometry::Camera camera(Math::makeVector(300.f, 0.f, 100.0f), Math::makeVector(0.f, 0.f, 0.f), 0.3f, ((double)visu.width()) / ((double)visu.height()), 1.0f);
		//camera.translateLocal(Math::makeVector(2000.f, 3000.f, 700.f));
		scene.setCamera(camera);
	}
	createGround(scene);
}




//4 000 000 de triangles ... en temps reel
void initEngine(Geometry::Scene & scene, Visualizer::Visualizer const& visu)
{
	Geometry::BoundingBox box(Math::makeVector(0.0f, 0.0f, 0.0f), Math::makeVector(0.0f, 0.0f, 0.0f));
	Geometry::Loader3ds loader(m_modelDirectory + "\\engine\\engine.3ds", m_modelDirectory + "\\engine\\Texture");
	// We remove the specular components of the materials...
	::std::vector<Geometry::Material*> materials = loader.getMaterials();
	for (auto it = materials.begin(), end = materials.end(); it != end; ++it)
	{
		//(*it)->setSpecular(RGBColor());
		(*it)->setAmbient((*it)->getAmbient() * 0.01);
		//(*it)->setEmissive(RGBColor());
	}

	for (size_t cpt = 0; cpt < loader.getMeshes().size(); ++cpt)
	{
		scene.add(*loader.getMeshes()[cpt]);
		box.update(*loader.getMeshes()[cpt]);
	}

	// 2.2 Adds point lights in the scene 
	Geometry::BoundingBox sb = scene.getBoundingBox();
	Math::Vector3f position = sb.max();
	Geometry::PointLight light1(position, RGBColor(1.0, 1.0, 1.0) * 1);
	scene.add(light1);
	position[1] = -position[1];
	Geometry::PointLight light2(position, RGBColor(1.0, 1.0, 1.0) * 1);
	scene.add(light2);
	{
		Geometry::Camera camera(Math::makeVector(30.f, 0.f, 10.0f), Math::makeVector(0.f, 0.f, 0.f), 0.3f, ((double)visu.width()) / ((double)visu.height()), 1.0f);
		//camera.translateLocal(Math::makeVector(2000.f, 3000.f, 700.f));
		scene.setCamera(camera);
	}
	createGround(scene);
}






////////////////////////////////////////////////////////////////////////////////////////////////////
/// \fn	void waitKeyPressed()
///
/// \brief	Waits until a key is pressed.
/// 		
/// \author	F. Lamarche, Université de Rennes 1
/// \date	03/12/2013
////////////////////////////////////////////////////////////////////////////////////////////////////
void waitKeyPressed()
{
	
	SDL_Event event;
	bool done = false;
	while (!done) {
		while (SDL_PollEvent(&event)) {
			switch (event.type) {
			case SDL_KEYDOWN:
				/*break;*/
			case SDL_QUIT:
				done = true;
				break;
			default:
				break;
			}
		}/*while*/
	}/*while(!done)*/
}




const std::string help_message =
"---REAL TIME RENDERING---\n"
"--Controls:\n"
"- Z Q S D Space LCtrl: Camera mouvement\n"
"- Arrows: Camera rotation\n"
"- LShift: Going faster\n"
"- CapsLock: Lock Going faster\n"
"- A: Launch global illumination (Stop with RETURN)\n"
"- I: Print infos (camera and visualizer)\n"
"- 0 - 9: set resolution scale\n"
"- F: Print FrameRate\n"
"- Rendering modes:\n"
"	- R: \"Realistic\" Phong rendering\n"
"	- B: Boxes of the acceleration structure\n"
"	- P: Depth\n"
"	- N: Normals\n"
"	- U: UV coordinates\n"
"	- K: Cartoon (TODO)\n"
"\n"
"- M: Mute (TODO)\n"
"- H: print Help\n"
;



void get_input(std::vector<SDL_Event> const& events,bool * keys, unsigned char & render_mode, int & scale, bool & rt)
{
	for(SDL_Event const& e: events)
	{
		if (e.type == SDL_KEYDOWN)
		{
			switch (e.key.keysym.sym)
			{
			case SDLK_UP:
				keys[0]=1;
				break;

			case SDLK_DOWN:
				keys[1] = 1;
				break;

			case SDLK_LEFT:
				keys[2] = 1;
				break;

			case SDLK_RIGHT:
				keys[3] = 1;
				break;
			case SDLK_z:
				keys[4] = 1;
				break;
			case SDLK_s:
				keys[5] = 1;
				break;
			case SDLK_d:
				keys[6] = 1;
				break;
			case SDLK_q:
				keys[7] = 1;
				break;
			case SDLK_SPACE:
				keys[8] = 1;
				break;
			case SDLK_LCTRL:
				keys[9] = 1;
				break;
			case SDLK_LSHIFT:
				keys[10] = 1;
				break;
			case SDLK_CAPSLOCK:
				keys[11] = ! keys[11];
				break;
			case SDLK_i:
				keys[12] = 1;
				break;
			case SDLK_f:
				keys[13] = !keys[13];
				break;
			case SDLK_b:
				render_mode = 0;
				break;
			case SDLK_r:
				render_mode = 1;
				break;
			case SDLK_u:
				render_mode = 2;
				break;
			case SDLK_n:
				render_mode = 3;
				break;
			case SDLK_p:
				render_mode = 4;
				break;
			case SDLK_a:
				rt = 0;
				break;
			case SDLK_h:
				std::cout << help_message << std::endl;
				break;
			default:
				break;
			}

			//NUMBERS
			if (e.key.keysym.sym >= SDLK_0 && e.key.keysym.sym <= SDLK_9)
			{
				scale = e.key.keysym.sym == SDLK_0 ? 10 : e.key.keysym.sym - SDLK_0;
			}
		}
		else if (e.type == SDL_KEYUP)
		{
			switch (e.key.keysym.sym)
			{
			case SDLK_UP:
				keys[0] = 0;
				break;

			case SDLK_DOWN:
				keys[1] = 0;
				break;

			case SDLK_LEFT:
				keys[2] = 0;
				break;

			case SDLK_RIGHT:
				keys[3] = 0;
				break;
			case SDLK_z:
				keys[4] = 0;
				break;
			case SDLK_s:
				keys[5] = 0;
				break;
			case SDLK_d:
				keys[6] = 0;
				break;
			case SDLK_q:
				keys[7] = 0;
				break;
			case SDLK_SPACE:
				keys[8] = 0;
				break;
			case SDLK_LCTRL:
				keys[9] = 0;
				break;
			case SDLK_LSHIFT:
				keys[10] = 0;
				break;
			default:
				break;
			}
		}
	}
}



////////////////////////////////////////////////////////////////////////////////////////////////////
/// \fn	int main(int argc, char ** argv)
///
/// \brief	Main entry-point for this application.
///
/// \author	F. Lamarche, Université de Rennes 1
/// \date	03/12/2013
///
/// \param	argc	Number of command-line arguments.
/// \param	argv	Array of command-line argument strings.
///
/// \return	Exit-code for the process - 0 for success, else an error code.
////////////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char ** argv)
{
	omp_set_num_threads(8);

	int scale = 5;

	// 1 - Initializes a window for rendering
	//Visualizer::Visualizer visu(2000, 2000, scale);// pour les ecrans 4K
	Visualizer::Visualizer visu(1000, 1000, scale);
	//Visualizer::Visualizer visu(1900, 1000, scale);
	//Visualizer::Visualizer visu(500, 500, scale);
	//Visualizer::Visualizer visu(300, 300, scale) ;
	//Visualizer::Visualizer visu(250, 250, scale) ;
	//Visualizer::Visualizer visu(200, 200, scale) ;
	//Visualizer::Visualizer visu(150, 150, scale) ;
	//Visualizer::Visualizer visu(100, 100, scale) ;

	// 2 - Initializes the scene
	//Geometry::Scene scene(&visu);
	//Geometry::bvh_scene scene(&visu);
	Geometry::turbo_bvh_scene scene(&visu);

	// 2.1 initializes the geometry (choose only one initialization)
	//init_test(scene, visu);
	//initDiffuse(scene, visu);
	//initDiffuse_dif(scene, visu);
	initDiffuse_surface_light(scene, visu);
	//initDiffuseSpecular(scene, visu) ;//custom
	//initSpecular(scene, visu) ;
	//initGuitar(scene, visu);
	//initDog(scene, visu);
	//initGarage(scene, visu);
	//initRobot(scene, visu);
	//initTemple(scene, visu);
	//initGraveStone(scene, visu);
	//initBoat(scene, visu);
	//initSombrero(scene, visu);
	//initTibetHouse(scene, visu);
	//initTibetHouseInside(scene, visu);
	//initMedievalCity(scene, visu);
	//initEngine(scene, visu);
	

	std::cout << "Building the acceleration structure" << std::endl;
	scene.pre_compute(8, 1, 1);
	std::cout << "Done!" << std::endl;
	

	////////////////////////////////////////////////////
	/// The shadow pre processing is currently desactivated,
	/// because it has some little ugly artifacts 
	////////////////////////////////////////////////////
	/*
	std::cout << "Pre computing the shadows" << std::endl;
	scene.pre_compute_shadows(10, 10, 10, 10);
	std::cout << "Done!" << std::endl;
	std::cout << scene.num_shadow << std::endl;
	*/

	// Shows stats
	scene.printStats();

	// 3 - Computes the scene
	unsigned int passPerPixel = 16;	// Number of rays per pixel 
	unsigned int subPixelSampling = 16;	// Antialiasing
										
	unsigned int maxBounce = 5;			// Maximum number of bounces

											//scene.setDiffuseSamples(16);
											//scene.setSpecularSamples(16);
	scene.setDiffuseSamples(1);
	scene.setSpecularSamples(1);
	//scene.setDiffuseSamples(32);
	//scene.setSpecularSamples(32);
	//scene.setDiffuseSamples(16);
	//scene.setSpecularSamples(16);
	//scene.setDiffuseSamples(4);
	//scene.setSpecularSamples(4);
	scene.setLightSamples(1);
	
	//real time option
	bool rt = 1;

	unsigned char render_mode = 1;


	std::cout << help_message << std::endl;


	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
	Geometry::Camera & cam = scene.m_camera;

	Math::Vector3f pos = cam.getPosition();
	Math::Vector3f dir = cam.getRay(0.5, 0.5).direction().normalized();

	//UP DOWN LEFT RIGHT Z S D Q SPACE LCTRL LSHIFT CAPSLOCK C F
	bool keys[14];
	for (int i = 0; i < 14; ++i)
	{
		keys[i] = 0;
	}
		


	const double speed = 2;
	const double angle_speed = 2;

	double forward(0), rightward(0), up(0), inclination, azimuth;

	Math::Vector3f speed_vec, cam_speed_vec;

	Math::Vector3f dir_sphere = spherical_coordinate(dir);

	inclination = dir_sphere[1];
	azimuth = dir_sphere[2];

		
	std::chrono::high_resolution_clock::time_point t2;
	unsigned int frame_count=0;
	while (1)
	{
		t2 = std::chrono::high_resolution_clock::now();
		get_input(visu.events, keys, render_mode, scale, rt);

		visu.update_scale(scale);

		if (rt)
		{

			std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);

			double dt = time_span.count();

			t1 = t2;


			if (keys[0])
			{
				inclination -= angle_speed * dt;
			}
			if (keys[1])
			{
				inclination += angle_speed * dt;
			}
			if (keys[2])
			{
				azimuth += angle_speed * dt;
			}
			if (keys[3])
			{
				azimuth -= angle_speed * dt;
			}
			if (keys[4])
			{
				forward += speed * dt;
			}
			if (keys[5])
			{
				forward -= speed * dt;
			}
			if (keys[6])
			{
				rightward += speed * dt;
			}
			if (keys[7])
			{
				rightward -= speed * dt;
			}
			if (keys[8])
			{
				up += speed * dt;
			}
			if (keys[9])
			{
				up -= speed * dt;
			}

			if (inclination > Math::pi)
			{

				inclination = Math::pi - 0.00000001;
			}
			else if (inclination < 0)
			{

				inclination = 0.00000001;
			}

			double sin_inc = sin(inclination);
			double cos_inc = cos(inclination);
			if (inclination > Math::pi)
			{
				cos_inc = -cos_inc;
			}
			dir = Math::makeVector(sin_inc * cos(azimuth), sin_inc * sin(azimuth), cos_inc);

			
			speed_vec = cam.m_right*rightward + cam.m_front*forward - cam.m_down*up;
			if (keys[10])
			{
				speed_vec = speed_vec * 15;
			}
			if (keys[11])
			{
				speed_vec = speed_vec * 15;
			}
			pos = pos + speed_vec;
			dir = dir + speed_vec;

			cam.update_both(pos + speed_vec, pos + dir);


			scene.realtime_compute(maxBounce, render_mode);


			forward = 0;
			rightward = 0;
			up = 0;

			//framerate

			if (keys[13] && !(frame_count % 50))
			{
				std::cout << "fr: " << 1.0 / dt << std::endl;
				//std::cout << "mean: " << mean_framerate << std::endl;
			}
			++frame_count;

			if (keys[12])
			{
				std::cout << "Camera position: " << cam.m_position << "\t Camera target: " << cam.m_target << std::endl;
				visu.print_info(std::cout);
				keys[12] = 0;
				
			}
		}
		else
		{
			std::cout<<
			"---All pass Mode---\n"
			"Press RETURN to stop\n"
			<< std::endl;
			scene.compute(maxBounce, render_mode, subPixelSampling, passPerPixel);
			rt = true;
		}
	}


	// 4 - waits until a key is pressed
	waitKeyPressed();
	
	return 0;
}
