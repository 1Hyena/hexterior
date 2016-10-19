#include <cmath>
#include <stdio.h>
#include <stdarg.h>
#include <array>
#include <vector>

#include "main.h"
#include "log.h"
#include "options.h"
#include "signals.h"
#include "allegro.h"
#include "camera.h"
#include "grid.h"

OPTIONS options;
SIGNALS signals;
ALLEGRO allegro;

std::vector<ALLEGRO_VERTEX> vertices;

int main(int argc, char **argv) {
    HEXAGON::size = 0.025;
    HEXAGON::gap  = 0.10;
    GRID grid;
    grid.fill(100);
    grid.generate_heightmap(2);
    
    size_t top_level = 7;
    for (size_t i=1777; i<3000; i+=1) {
        size_t peel, index;
        HEXAGON::vortex_to_polar(i, &peel, &index);
        //if (grid[peel][index].is_master(top_level)) continue;
        grid[peel][index].highlight = true;

        std::map<size_t, size_t> masters;
        grid[peel][index].to_masters(top_level, &masters);
        //printf("Checking %lu masters of %lu:%lu (%d, %d), vortex %lu.\n", masters.size(), peel, index, grid[peel][index].get_x(), grid[peel][index].get_y(), grid[peel][index].get_vortex());
        bool break_all = false;
        for (const auto & m : masters) {
            HEXAGON::vortex_to_polar(m.first, &peel, &index);
            if (grid[peel][index].is_master(top_level) == false) {
                printf("vortex %lu (%lu:%lu) is not master!\n", grid[peel][index].get_vortex(), peel, index);
                break_all = true;
            } else grid[peel][index].highlight = true;
        }
        if (break_all) break;
        break;
    }
    
    /*
    std::map<size_t, size_t> masters;
    grid[2][0].to_masters(1, &masters);
    for (const auto & m : masters) {
        size_t peel, index;
        HEXAGON::vortex_to_polar(m.first, &peel, &index);
        grid[peel][index].highlight = true;
    }
    */
    grid.to_vertices(&vertices);
    //printf("%lu\n", vertices.size());

    //grid.populate(0);
    //HEXAGON hex;
    //printf("%f \n", grid[0][10].get_size());//"%f\n", grid[0][1].get_size());
    if (!init(argc, argv)) {
        deinit(true);
        fprintf(stderr, "%s: failed to initialize.\n", argv[0]);
        return EXIT_FAILURE;
    }
    if (options.exit_flag) return deinit(true);

    log("\x1B[1;32mHexterior is ready to rock!\x1B[0m\x1B]0;Hexterior\a");

    CAMERA camera;

    int redraw = 0;
    al_start_timer(allegro.timer);
    while (true) {
        ALLEGRO_EVENT event;

        al_wait_for_event(allegro.queue, &event);
        if (event.type == ALLEGRO_EVENT_DISPLAY_CLOSE) break;
        else if (event.type == ALLEGRO_EVENT_DISPLAY_RESIZE) {
            al_acknowledge_resize(allegro.display);
        }
        else if (event.type == ALLEGRO_EVENT_KEY_DOWN) {
            if (event.keyboard.keycode == ALLEGRO_KEY_ESCAPE) break;
            if (event.keyboard.keycode == ALLEGRO_KEY_SPACE) {
                allegro.controls++;
                allegro.controls %= 3;
            }
            allegro.key[event.keyboard.keycode] = 1;
            allegro.keystate[event.keyboard.keycode] = 1;
        }
        else if (event.type == ALLEGRO_EVENT_KEY_UP) {
            // In case a key gets pressed and immediately released, we will still
            // have set ex.key so it is not lost.
            allegro.keystate[event.keyboard.keycode] = 0;
        }
        else if (event.type == ALLEGRO_EVENT_TIMER) {
            bool terminate = false;
            // All non-fatal signals must be blocked before the resulting atomic
            // flags are checked because otherwise the signal could go lost if
            // it appears after atomic flag checking and before sockets.wait.
            signals.block();
            while (int sig = signals.next()) {
                char *str = strsignal(sig);
                switch (sig) {
                    case SIGINT :
                    case SIGTERM:
                    case SIGQUIT: terminate = true;
                    default     : log("Caught signal %d (%s).", sig, str ? str : "NULL"); break;
                }
            }
            signals.unblock();
            if (terminate) break;

            camera.handle_input();
            redraw = 1;

            // Reset keyboard state for keys not held down anymore.
            for (int i = 0; i < ALLEGRO_KEY_MAX; i++) {
                if (allegro.keystate[i] == 0) allegro.key[i] = 0;
            }
            allegro.mouse_dx = 0;
            allegro.mouse_dy = 0;
        }
        else if (event.type == ALLEGRO_EVENT_MOUSE_BUTTON_DOWN) {
            allegro.button[event.mouse.button] = 1;
        }
        else if (event.type == ALLEGRO_EVENT_MOUSE_BUTTON_UP) {
            allegro.button[event.mouse.button] = 0;
        }
        else if (event.type == ALLEGRO_EVENT_MOUSE_AXES) {
            allegro.mouse_dx += event.mouse.dx;
            allegro.mouse_dy += event.mouse.dy;
        }

        if (redraw  && al_is_event_queue_empty(allegro.queue)) {
            int fps;
            if ( (fps = allegro.calculate_fps()) == -1) continue;

            al_clear_to_color(al_map_rgb(0,0,0));
            camera.draw_scene();
            al_draw_textf(allegro.font, al_map_rgb(0,255,0), 0, 0, 0, "FPS: %3d", fps);
            al_flip_display();
            redraw = 0;
        }
    }

    return deinit(false);
}

bool init(int argc, char **argv) {
    std::string version;
    version.append("Hexterior version ");
    version.append(HEXTERIOR_VERSION);
    version.append(" Copyright (C) 2016 Erich Erstu");

    if (!options.init(&log_options, argc, argv, version.c_str())
    ||  !signals.init(&log_signals)) return false;

    if (options.exit_flag) return true;

    allegro.preset_max_fps(60.0);
    allegro.preset_resolution(640, 360);
    if (!allegro.init(&log_allegro)) return false;
/*
    size_t w = heightmap.size();
    for (size_t x=0; x<w; ++x) {
        size_t h = heightmap[x].size();
        for (size_t y=0; y<h; ++y) {
            heightmap[x][y].init(x, y, w, h);
        }
    }
*/
    return true;
}

int deinit(bool silent) {
    vertices.clear();

    int result = EXIT_SUCCESS;
    if (!allegro.deinit()) result = EXIT_FAILURE;

    if (!silent) {
        if (result == EXIT_FAILURE)  log("\x1B[1;31mHexterior closes with errors!\x1B[0m");
        else if (!options.exit_flag) log("\x1B[1;33mNormal termination of Hexterior.\x1B[0m");
    }

    return result;
}

void vlog(const char *format, ...) {
    if (!options.verbose) return;
    char buffer[256];
    va_list args;
    va_start(args, format);
    vsnprintf(buffer, sizeof(buffer), format, args);
    va_end(args);
    log_text(buffer);
}

