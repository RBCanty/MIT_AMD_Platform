import mcn_logging


# Create and bind logger for the AMDES's systems
system_log = mcn_logging.make_logger("DEBUG", "MCN Logger")
system_log.addHandler(mcn_logging.make_console_handler())
log_file = mcn_logging.make_logfile_path("mcn_logfile.log")
mcn_logging.bind_handler(system_log, log_file)

# Attach methods to system_log object such that the metric monitor
#   can make adjustments for the system being on or off
UP = "__SYSTEM__::MCN_OPERATION IS *ON*"
DOWN = "__SYSTEM__::MCN_OPERATION IS *OFF*"
system_log.mark_up = lambda: system_log.info(UP)
system_log.mark_down = lambda: system_log.info(DOWN)


if __name__ == '__main__':
    print(log_file)
    system_log.info("Hello World")