import os
from configparser import ConfigParser, ExtendedInterpolation

configsPath = os.path.realpath(__file__).rstrip("/__init__.py")

pConfig = ConfigParser(interpolation = ExtendedInterpolation())
pConfig.read(os.path.join(configsPath, "processes.ini"))

dConfig = ConfigParser(interpolation = ExtendedInterpolation())
dConfig.read(os.path.join(configsPath, "directories.ini"))

loggingConfigPath = lConfigPath = LOGGING_CONFIG = os.path.join(configsPath, "logging.ini")

sConfig = ConfigParser(interpolation = ExtendedInterpolation())
sConfig.read(os.path.join(configsPath, "setup.ini"))
