import mcn_logging


# Create and bind logger for the metric monitoring
metric_log = mcn_logging.make_logger("DEBUG", "Metric Logger")
metric_log.addHandler(mcn_logging.make_console_handler())
metrics_log_file = mcn_logging.make_logfile_path("analytics.log")
mcn_logging.bind_handler(metric_log, metrics_log_file)


if __name__ == '__main__':
    print(metrics_log_file)
    metric_log.info("Hello World")