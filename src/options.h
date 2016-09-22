#include <string>
#include <getopt.h>
#include <sys/types.h>
#include <dirent.h>

class OPTIONS {
    public:
    OPTIONS()  {}
    ~OPTIONS() {}

    int verbose             =  0;
    int exit_flag           =  0;
    std::string name        = "";
    std::string datadir     = ".";

    inline void print_usage (FILE* stream) {
        fprintf (stream, "Usage: %s [options]\n", name.c_str());
        fprintf (stream, "Options:\n");
        fprintf (stream,
            "      --brief           Print brief messages (default).\n"
            "  -d  --datadir         Specify data directory.\n"
            "  -h  --help            Display this usage information.\n"
            "      --verbose         Print verbose messages.\n"
            "  -v  --version         Show version information.\n"
        );
    }

    inline bool init(void (*log_function)(const char *p_fmt, ...), int argc, char **argv, const char *version) {
        if (log_function) log = log_function;

        int c;
        name = argv[0];
        while (1) {
            static struct option long_options[] = {
                // These options set a flag:
                {"brief",               no_argument,         &verbose,           0 },
                {"verbose",             no_argument,         &verbose,           1 },
                // These options don't set a flag. We distinguish them by their indices:
                {"help",                no_argument,              0,            'h'},
                {"version",             no_argument,              0,            'v'},
                {"datadir",             required_argument,        0,            'd'},
                {0,                     0,                        0,             0 }
            };

            // getopt_long stores the option index here.
            int option_index = 0;

            c = getopt_long(argc, argv, "d:hv", long_options, &option_index);

            /* Detect the end of the options. */
            if (c == -1) break;

            switch (c) {
                case 0:
                    {
                        // If this option sets a flag do nothing else.
                        if (long_options[option_index].flag != 0) break;
                        std::string buf="option ";
                        buf.append(long_options[option_index].name);
                        if (optarg) {
                            buf.append(" with arg ");
                            buf.append(optarg);
                        }
                        log(buf.c_str());
                        break;
                    }
                case 'd':
                    datadir = optarg;
                    break;
                case 'h': print_usage(stdout); exit_flag = 1; break;
                case 'v':
                    printf("%s\n", version);
                    exit_flag = 1;
                    break;
                case '?':
                    // getopt_long already printed an error message.
                    break;
                default: return false;
            }
        }

        if (exit_flag) return true;

        if (datadir.length() > 0) {
            DIR* dir = opendir(datadir.c_str());
            if (dir) closedir(dir);
            else if (ENOENT == errno) {
                log("Data directory %s does not exist.", datadir.c_str());
                return false;
            }
            else {
                log("Could not open data directory %s.", datadir.c_str());
                return false;
            }
        }
        else {
            log("Data directory must be defined.");
            return false;
        }

        while (optind < argc) log("unidentified argument: %s", argv[optind++]);
        return true;
    }

    private:
    inline static void drop_log(const char *, ...) {}

    void (*log)(const char *p_fmt, ...) = drop_log;
};

