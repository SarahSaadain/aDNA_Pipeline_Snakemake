import yaml

from common.common_config_enumerations import ConfigSettings
from types import SimpleNamespace
from enum import Enum

class ConfigManager:
    def __init__(self, config_dict):
        self.config = config_dict

    def _dict_to_namespace(self, d) -> SimpleNamespace:
        if isinstance(d, dict):
            return SimpleNamespace(**{k: self._dict_to_namespace(v) for k, v in d.items()})
        elif isinstance(d, list):
            return [self._dict_to_namespace(i) for i in d]
        else:
            return d

    def get_value(self, *path, default=None) -> any:
        """
        Generic deep get for arbitrary nested keys, returning dict or default.
        """
        node = self.config
        for key in path:
            if not isinstance(node, dict):
                return default
            node = node.get(key, default)
            if node == default:
                return default
        return node

    def get_step_setting(self, *path) -> SimpleNamespace:
        """
        Get settings under the pipeline root with attribute access.

        Usage:
            get_step_setting("run_raw_reads_processing", "adapter_removal")
        """
        node = self.config.get(ConfigSettings.PIPELINE_CONFIG.value, {})
        for key in path:
            node = node.get(key, {})
        return self._dict_to_namespace(node)
    
    def get_species_list(self) -> list:
        """
        Get a list of species from the config.
        """
        species = self.config.get("species", {})
        return list(species.keys())
    
    def get_log_level(self) -> str:
        """
        Get the logging level from the config.
        """
        return self.get_value(ConfigSettings.PIPELINE_CONFIG.value, "log_level", default="INFO").upper()

    def get_species_setting(self, species_key) -> SimpleNamespace:
        """
        Get species configuration as a SimpleNamespace.
        """
        species = self.config.get("species", {})
        return self._dict_to_namespace(species.get(species_key, {}))

    def get_compare_species_setting(self, compare_key) -> SimpleNamespace:
        """
        Get compare_species configuration as a SimpleNamespace.
        """
        compares = self.config.get("compare_species", {})
        return self._dict_to_namespace(compares.get(compare_key, {}))
    
    def get_threads(self) -> int:
        """
        Get the number of threads from the config.
        """
        return self.get_value(ConfigSettings.PIPELINE_CONFIG.value, "threads", default=5)
    
    def get_project_path(self):
        """
        Get the project path from the config.
        """
        return self.config.get("path_adna_project", "")
    
    def get_pipeline_settings(self, stage_key: Enum, process_key: Enum = None, species: str = None) -> dict:
        
        stage_key_str = stage_key.value

         # Initialize species settings dictionary
        species_settings = {}
    
        # Construct path for general settings (under 'processing')
        if process_key is None:
            # If process_key is not provided, use stage_key to construct the path
            general_path = [ConfigSettings.PIPELINE_CONFIG.value, stage_key_str]
        else:

            process_key_str = process_key.value

            # If process_key is provided, use it to construct the path
            general_path = [ConfigSettings.PIPELINE_CONFIG.value, stage_key_str, process_key_str]

            # If species is provided, construct path and retrieve species-specific settings
            if species:
                species_path = [ConfigSettings.SPECIES_CONFIG.value, species, ConfigSettings.PIPELINE_CONFIG.value, stage_key_str, process_key_str]
                species_settings = self.get_value(*species_path, default={})

        # Retrieve general settings
        general_settings = self.get_value(*general_path, default={})

        # Merge settings: general_settings first, then species_settings to overwrite
        combined_settings = self._deep_merge(general_settings, species_settings)

        # Ensure the result is a dictionary
        return combined_settings if isinstance(combined_settings, dict) else {}

    def is_process_stage_enabled(self, stage_key: Enum, default: bool = True) -> bool:
        """
        Check if a specific process step is enabled for a given species.
        
        Args:
            process_key (Enum): The process key to check.
            species (str): The species to check the process for. If None, checks general settings.
        
        Returns:
            bool: True if the process step is enabled, False otherwise.
        """
        settings = self.get_pipeline_settings(stage_key)
        return settings.get("execute", default)
    
    def is_process_step_enabled(self, stage_key: Enum, step_key: Enum, species: str = None, default: bool = True) -> bool:
        """
        Check if a specific subprocess is enabled for a given process step.
        
        Args:
            process_key (Enum): The process key to check.
            subprocess_key (str): The subprocess key to check.
        
        Returns:
            bool: True if the subprocess is enabled, False otherwise.
        """
        settings = self.get_pipeline_settings(stage_key, step_key, species)
        return settings.get("execute", default)
    
    def _deep_merge(self, dict1, dict2):
        result = dict1.copy()
        for key, value in dict2.items():
            if (
                key in result 
                and isinstance(result[key], dict) 
                and isinstance(value, dict)
            ):
                result[key] = self._deep_merge(result[key], value)
            else:
                result[key] = value
        return result


_config = None  # Private module-level variable

def load_config(config_file) -> dict:
    """Loads the configuration from the specified YAML file."""
    global _config
    
    with open(config_file, 'r') as f:
        _config = yaml.safe_load(f)
        
    return _config