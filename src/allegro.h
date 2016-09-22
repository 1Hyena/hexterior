#ifndef ALLEGRO_H_15_09_2016
#define ALLEGRO_H_15_09_2016

#include <allegro5/allegro.h>
#include <allegro5/allegro_image.h>
#include <allegro5/allegro_primitives.h>
#include <allegro5/allegro_color.h>
#include <allegro5/allegro_font.h>

class ALLEGRO {
    public:
    ALLEGRO()  {}
    ~ALLEGRO() {}

    ALLEGRO_DISPLAY *display   = nullptr;
    ALLEGRO_TIMER *timer       = nullptr;
    ALLEGRO_EVENT_QUEUE *queue = nullptr;
    ALLEGRO_FONT *font         = nullptr;

    // controls sensitivity
    double mouse_look_speed = 0.03;
    double movement_speed = 0.05;

    // keyboard and mouse state
    int button[10];
    int key[ALLEGRO_KEY_MAX];
    int keystate[ALLEGRO_KEY_MAX];
    int mouse_dx, mouse_dy;
    int controls=0;


    bool init(void (*log_function)(const char *p_fmt, ...)) {
        if (log_function) log = log_function;

        if (calc_fps_max_fps <= 0.0) {
            log("FPS must be a positive number.");
            return false;
        }

        if (screen_w < 1 || screen_h < 1) {
            log("Invalid display resolution.");
            return false;
        }

        if (!al_install_system(ALLEGRO_VERSION_INT, nullptr)) {
            log("Could not install system.");
            return false;
        }
        sys_installed = true;

        if (!al_init_font_addon()) {
            log("Could not initialize font addon.");
            return false;
        }

        if (!al_init_primitives_addon()) {
            log("Could not initialize primitives addon.");
            return false;
        }
        
        if (!al_init_image_addon()) {
            log("Could not initialize image addon.");
            return false;
        }

        if (!al_install_keyboard()) {
            log("Could not install keyboard.");
            return false;
        }

        if (!al_install_mouse()) {
            log("Could not install mouse.");
            return false;
        }

        al_set_new_display_option(ALLEGRO_SAMPLE_BUFFERS, 1, ALLEGRO_SUGGEST);
        al_set_new_display_option(ALLEGRO_SAMPLES, 8, ALLEGRO_SUGGEST);
        al_set_new_display_option(ALLEGRO_DEPTH_SIZE, 16, ALLEGRO_SUGGEST);
        al_set_new_display_flags(ALLEGRO_RESIZABLE);
        if ((display = al_create_display(screen_w, screen_h)) == nullptr) {
            log("Error creating display.");
            return false;
        }

        if ((font = al_create_builtin_font()) == nullptr) {
            log("Failed to create a builtin font.");
            return false;
        }

        timer = al_create_timer(1.0 / calc_fps_max_fps);
        queue = al_create_event_queue();
        al_register_event_source(queue, al_get_keyboard_event_source());
        al_register_event_source(queue, al_get_mouse_event_source());
        al_register_event_source(queue, al_get_display_event_source(display));
        al_register_event_source(queue, al_get_timer_event_source(timer));

        return true;
    }

    bool deinit() {
        if (al_is_system_installed()) {
            if (queue)   al_destroy_event_queue(queue);
            if (timer)   al_destroy_timer(timer);
            if (display) al_destroy_display(display);
            if (font)    al_destroy_font(font);

            al_uninstall_system();
        }
        else if (sys_installed) log("Warning! System is already uninstalled.");

        return true;
    }

    int calculate_fps() {
        if (calc_fps_first) {
            calc_fps_first = false;
            calc_fps_old_time = al_get_time();
            return -1;
        }

        int rec_times    = 0;
        int max_times    = round_int(calc_fps_max_fps);
        double new_time  = al_get_time();
        double delta     = new_time - calc_fps_old_time;
        calc_fps_delta_sum += delta;
        calc_fps_old_time   = new_time;
        double p  = calc_fps_delta_sum * max_times;
        rec_times = round_int(p);

        if (calc_fps_times > rec_times) {
            return -1;
        }
        calc_fps_times++;

        int fps = 0;
        if (calc_fps_delta_sum >= 1.0 || calc_fps_times>=max_times) {
            fps = calc_fps_times;
            calc_fps_old_fps = fps;
            calc_fps_times=0;
            calc_fps_delta_sum=0.0;
        }
        else {
            if (calc_fps_old_fps == -1) fps = calc_fps_times;
            else                        fps = calc_fps_old_fps;
        }

        return fps;
    }

    // Set up a perspective transform. We make the screen span
    // 2 vertical units (-1 to +1) with square pixel aspect and the camera's
    // vertical field of view. Clip distance is always set to 1.
    void setup_3d_projection() {
        ALLEGRO_TRANSFORM projection;
        double dw = al_get_display_width(display);
        double dh = al_get_display_height(display);
        double f;
        double vertical_field_of_view = 60.0 * ALLEGRO_PI / 180.0;
        al_identity_transform(&projection);
        al_translate_transform_3d(&projection, 0, 0, -1);
        f = tan(vertical_field_of_view / 2);
        al_perspective_transform(&projection, -1 * dw / dh * f, f, 1, f * dw / dh, -f, 1000);
        al_use_projection_transform(&projection);
    }

    void preset_max_fps   (double max_fps)                 { if (sys_installed) return; calc_fps_max_fps = max_fps; }
    void preset_resolution(unsigned int w, unsigned int h) { if (sys_installed) return; screen_w = w; screen_h = h; }

    private:
    bool sys_installed = false;

    int    calc_fps_times     = 0;
    double calc_fps_old_time  = 0.0;
    double calc_fps_delta_sum = 0.0;
    int    calc_fps_old_fps   = -1;
    bool   calc_fps_first     = true;
    double calc_fps_max_fps   = 60.0;
    int    screen_w           = 640;
    int    screen_h           = 480;

    static void drop_log(const char *, ...) {}

    void (*log)(const char *p_fmt, ...) = drop_log;

    int round_int(double r) {return (r > 0.0) ? (r + 0.5) : (r - 0.5);}
};

#endif

