#include "allegro.h"
#include "vector_3d.h"

extern ALLEGRO allegro;
extern std::vector<ALLEGRO_VERTEX> vertices;

class CAMERA {
    public:
    CAMERA() {
        xaxis.x = 1;
        yaxis.y = 1;
        zaxis.z = 1;
        position.y = 2;
        vertical_field_of_view = 60.0 * ALLEGRO_PI / 180.0;
    };

   ~CAMERA() {};

    VECTOR_3D position;
    VECTOR_3D xaxis; // This represent the direction looking to the right.
    VECTOR_3D yaxis; // This is the up direction.
    VECTOR_3D zaxis; // This is the direction towards the viewer ('backwards').
    double vertical_field_of_view; // In radians.

    void draw_scene();
    void handle_input();
    void move_along_ground(double right, double forward);
    void rotate_around_axis(VECTOR_3D axis, double radians);
    VECTOR_3D get_ground_forward_vector();
    VECTOR_3D get_ground_right_vector();
};

inline void CAMERA::handle_input() {
    double x = 0, y = 0;
    double xy;
    if (allegro.key[ALLEGRO_KEY_A] || allegro.key[ALLEGRO_KEY_LEFT])  x = -1;
    if (allegro.key[ALLEGRO_KEY_S] || allegro.key[ALLEGRO_KEY_DOWN])  y = -1;
    if (allegro.key[ALLEGRO_KEY_D] || allegro.key[ALLEGRO_KEY_RIGHT]) x = 1;
    if (allegro.key[ALLEGRO_KEY_W] || allegro.key[ALLEGRO_KEY_UP])    y = 1;

    /* Change field of view with Z/X. */
    if (allegro.key[ALLEGRO_KEY_Z]) {
        double m = 20 * ALLEGRO_PI / 180;
        vertical_field_of_view -= 0.01;
        if (vertical_field_of_view < m) vertical_field_of_view = m;
    }
    if (allegro.key[ALLEGRO_KEY_X]) {
        double m = 120 * ALLEGRO_PI / 180;
        vertical_field_of_view += 0.01;
        if (vertical_field_of_view > m) vertical_field_of_view = m;
    }

    /* In FPS style, always move the camera to height 2. */
    if (allegro.controls == 0) {
        if (position.y > 2) position.y -= 0.1;
        if (position.y < 2) position.y  = 2;
    }

    /* Set the roll (leaning) angle to 0 if not in airplane style. */
    /*if (allegro.controls == 0 || allegro.controls == 2) {
        double roll = get_roll(&allegro.camera);
        camera_rotate_around_axis(&allegro.camera, allegro.camera.zaxis, roll / 60);
    }*/

    /* Move the camera, either freely or along the ground. */
    xy = sqrt(x * x + y * y);
    if (xy > 0) {
        x /= xy;
        y /= xy;
        if (allegro.controls == 0) {
            move_along_ground(allegro.movement_speed * x, allegro.movement_speed * y);
        }
        /*if (allegro.controls == 1 || allegro.controls == 2) {
            camera_move_along_direction(&allegro.camera, allegro.movement_speed * x,
            allegro.movement_speed * y);
        }*/
    }

    /* Rotate the camera, either freely or around world up only. */
    if (allegro.button[1]) {
        if (allegro.controls == 0 || allegro.controls == 2) {
        VECTOR_3D up(0, 1, 0);
            rotate_around_axis(xaxis, -allegro.mouse_look_speed * allegro.mouse_dy);
            rotate_around_axis(up,    -allegro.mouse_look_speed * allegro.mouse_dx);
        }
        /*if (allegro.controls == 1) {
        camera_rotate_around_axis(&allegro.camera, allegro.camera.xaxis,
        -allegro.mouse_look_speed * allegro.mouse_dy);
        camera_rotate_around_axis(&allegro.camera, allegro.camera.zaxis,
        -allegro.mouse_look_speed * allegro.mouse_dx);
        }*/
    }
}

inline void CAMERA::draw_scene() {
    /* We save Allegro's projection so we can restore it for drawing text. */
    ALLEGRO_TRANSFORM projection = *al_get_current_projection_transform();
    ALLEGRO_TRANSFORM t;

    allegro.setup_3d_projection();

    /* We use a depth buffer. */
    al_set_render_state(ALLEGRO_DEPTH_TEST, 1);
    al_clear_depth_buffer(1);

    /* Construct a transform corresponding to our camera. This is an inverse
    * translation by the camera position, followed by an inverse rotation
    * from the camera orientation.
    *//*
    float x = 0.0;
    float y = 2.0;
    float z = 0.0; 
    //float x_x = 1.0, x_y = 0.0, x_z = 0.0;
    float y_x = 0.0, y_y = 1.0, y_z = 0.0;
    float z_x = 0.0, z_y = 0.0, z_z = 1.0;
    al_build_camera_transform(&t, x, y, z, x - z_x, y - z_y, z - z_z, y_x, y_y, y_z);*/
    al_build_camera_transform(&t,
        position.x,
        position.y,
        position.z,
        position.x - zaxis.x,
        position.y - zaxis.y,
        position.z - zaxis.z, 
        yaxis.x, yaxis.y, yaxis.z
    );
    al_use_transform(&t);
/*
    ALLEGRO_VERTEX vs[3];
    vs[0].x =   0; vs[0].y = -0.2; vs[0].z = 12; vs[0].color = al_map_rgb_f(0, 1, 1);
    vs[1].x = -12; vs[1].y = -0.2; vs[1].z = -8; vs[1].color = al_map_rgb_f(1, 0, 1);
    vs[2].x =  12; vs[2].y = -0.2; vs[2].z = -8; vs[2].color = al_map_rgb_f(1, 1, 0);*/
    /*
    vs[3].x =   0; vs[3].y = -0.2; vs[3].z =  9; vs[3].color = al_map_rgb_f(1, 0, 0);
    vs[4].x =  -9; vs[4].y = -0.2; vs[4].z = -6; vs[4].color = al_map_rgb_f(0, 1, 0);
    vs[5].x =   9; vs[5].y = -0.2; vs[5].z = -6; vs[5].color = al_map_rgb_f(0, 0, 1);*/
    //al_draw_prim(vs, nullptr, nullptr, 0, 3, ALLEGRO_PRIM_TRIANGLE_LIST);

    al_draw_prim(&vertices[0], nullptr, nullptr, 0, vertices.size(), /*ALLEGRO_PRIM_POINT_LIST*/ALLEGRO_PRIM_TRIANGLE_LIST);

    /* Restore projection. */
    al_identity_transform(&t);
    al_use_transform(&t);
    al_use_projection_transform(&projection);
    al_set_render_state(ALLEGRO_DEPTH_TEST, 0);
}

/* Like camera_move_along_direction but moves the camera along the ground plane
 * only.
 */
inline void CAMERA::move_along_ground(double right, double forward) {
    VECTOR_3D f = get_ground_forward_vector();
    VECTOR_3D r = get_ground_right_vector();
    position.x += f.x * forward + r.x * right;
    position.z += f.z * forward + r.z * right;
}

/* Rotate the camera around the given axis. */
inline void CAMERA::rotate_around_axis(VECTOR_3D axis, double radians) {
    ALLEGRO_TRANSFORM t;
    al_identity_transform(&t);
    al_rotate_transform_3d(&t, axis.x, axis.y, axis.z, radians);
    al_transform_coordinates_3d(&t, &yaxis.x, &yaxis.y, &yaxis.z);
    al_transform_coordinates_3d(&t, &zaxis.x, &zaxis.y, &zaxis.z);

    /* Make sure the axes remain orthogonal to each other. */
    zaxis.normalize();
    xaxis = VECTOR_3D::cross_product(yaxis, zaxis);
    xaxis.normalize();
    yaxis = VECTOR_3D::cross_product(zaxis, xaxis);
}

/* Get a vector with y = 0 looking in the opposite direction as the camera z
 * axis. If looking straight up or down returns a 0 vector instead.
 */
inline VECTOR_3D CAMERA::get_ground_forward_vector() {
    VECTOR_3D move(zaxis);
    move.mul(-1);
    move.y = 0;
    move.normalize();
    return move;
}

/* Get a vector with y = 0 looking in the same direction as the camera x axis.
 * If looking straight up or down returns a 0 vector instead.
 */
inline VECTOR_3D CAMERA::get_ground_right_vector() {
    VECTOR_3D move(xaxis);
    move.y = 0;
    move.normalize();
    return move;
}

