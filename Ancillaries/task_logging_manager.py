import mcn_logging


# Create and bind a logger for the task manager (to keep clutter off of the main system logger)
tasking_log = mcn_logging.make_logger("DEBUG", "Tasking Logger")
tasking_log_file = mcn_logging.make_logfile_path("task_manager.log")
mcn_logging.bind_handler(tasking_log, tasking_log_file)


if __name__ == '__main__':
    print(tasking_log_file)
    tasking_log.info("Hello World")
