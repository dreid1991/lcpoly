// Povray scene description
#include "colors.inc"
#include "transforms.inc"
#include "textures.inc"
#include "metals.inc"

// camera and lights:
#declare H = 0; // height of camera location and its focus point
#declare lookat_location = <0,0,H>;
#declare camera_location = <0,40,H>;
#declare camera_rotation = <0,0,-180>;
#declare light1_position = <0,40,0>;
#declare light2_position = <0,0,40>;
#declare length = 25; // image width

#declare backgroundcolor = <1,1,1>;
// for object texture:
#declare transmit_l = 0.5;
#declare ambient_l = 0.3;
#declare diffuse_l = 0.9;
#declare brilliance_l = 1.0;
#declare phong_l = 0.5;

#declare White = rgb <0,1,0.4>;

//############ settings for camera/light and brightness ###################
global_settings
{
    assumed_gamma 1.8  // for device brightness correction
}
//-------------------------------

camera
{
    orthographic  // projection type: perspective(default) | orhographic | spherical
    location    camera_location

    right   -x*length                           // horizontal vector in image plane 
    up       y*length*image_height/image_width  // vertical vector in image plane scaled by 
                                                //   image height and width defined in RES.ini
    sky     <0,0,1>     // point z-axis up 
    rotate      camera_rotation
    look_at     lookat_location
}

// light source 1:
light_source
{
    light1_position     // light's position
    color rgb <1,1,1>   // light's color
    rotate      camera_rotation // rotate light with camera
    parallel    // light type modifiers: spotlight | shadowless | cylider | parallel
}

// light source 2:
light_source
{
    light2_position
    color rgb <1,1,1>
//    point_at  <0,0,0>
}

// background:
background { backgroundcolor }


//##################### macros definitions ################################
// define macro for object texture:
#macro object_texture( Color )
    texture 
    {
//        NewColor=Color*Color.blue
        pigment 
        {
            color Color
//            color Color/(-(Color-1)*(Color-1)+1.1)
//            color Color/(-(Color-1)*(Color-1)*(Color-1)*(Color-1)+1.1)
            transmit transmit_l
        }
        finish 
        {
//            ambient ambient_l
            diffuse diffuse_l
//            brilliance brilliance_l
            phong phong_l
        }
    }
#end

// define an ellipsoidal partical:
#macro ellipsoid( Position, Rotation, Color, Scale )
    sphere
    {
        <0,0,0>       // generated in <x,y,z>
        0.5           // radius of sphere
        scale Scale   // scale sphere by <Sx,Sy,Sz>
        // rotate ellipsoid by alinging Axis1 to Axis2:
        Reorient_Trans(<1,0,0>, Rotation)   // part of "transforms.inc"
        translate Position    // "translate <x,y,z>" will move an ellipsoid from current position <x0,y0,z0> to new one <x0+x,y0+y,z0+z>
        object_texture(Color) // specify texture using previously defined macro
    }
#end

// define a spherocylindrical particle:
#macro spherocylinder( Position, Rotation, Color, Radius, Length )
    merge
    {
      // if z-axis is longest:
//        cylinder{ <0,0,-0.5*Length>, <0,0,0.5*Length>, Radius }
//        sphere  { <0,0,-0.5*Length> Radius }
//        sphere  { <0,0, 0.5*Length> Radius }
      // if x-axis is longest:
        cylinder{ <-0.5*Length,0,0>, <0.5*Length,0,0>, Radius }
        sphere  { <-0.5*Length,0,0> Radius }
        sphere  { <0.5*Length,0,0> Radius }
        Reorient_Trans(<1,0,0>, Rotation)
        translate Position
        object_texture(Color)
    }
#end

// define a spheroplatelet particle:
#macro spheroplatelet( Position, Rotation, Color, Major, Minor )
    merge
    {
        torus   // assume x- and y-axes are long
        { Major, Minor
          rotate  <90,0,0> }

        cylinder{ <0,0,-Minor>, <0,0,Minor>, Major }

        Reorient_Trans(<1,0,0>, Rotation)
        translate Position
        object_texture(Color)
    }
#end


//################## particle list ########################################
