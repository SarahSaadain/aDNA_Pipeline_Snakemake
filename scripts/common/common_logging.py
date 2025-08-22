import logging

#####################
# Log
#####################

# --- Logging Setup (EARLY) ---
LOG_FORMAT = '[%(asctime)s] [%(levelname)s] %(message)s'
LOG_DATE_FORMAT = '%Y-%m-%d %H:%M:%S (%Z)'

logging.basicConfig(  # Basic config ASAP (for fallback)
    level=logging.INFO,
    format=LOG_FORMAT,
    datefmt=LOG_DATE_FORMAT,
    handlers=[logging.StreamHandler()]  # Only console for now
)

def print_command(subprocess_command: list):  # 🚀 Used for subprocess command execution
    command = ' '.join(subprocess_command)
    logging.info(f"🚀  {command}")

def print_info(message: str):
    logging.info(f"ℹ️  {message}")

def print_error(message: str):
    logging.error(f"❌ {message}")

def print_success(message: str):
    logging.info(f"✅ {message}")

def print_skipping(message: str):
    logging.info(f"⏩ Skipping: {message}")

def print_warning(message: str):
    logging.warning(f"⚠️ {message}")

def print_debug(message: str):
    logging.debug(f"🐞 {message}")

import logging

def print_step_execution(message: str):
    emoji = "🛠️"
    formatted_message = f"{emoji}  {message.upper()}  {emoji}"
    separator = "─" * (len(formatted_message) + 2)
    
    logging.info("")  # Optional line break before
    logging.info(separator)
    logging.info(f"│ {formatted_message} │")
    logging.info(separator)

def print_species_execution(message: str):
    
    # animal
    emoji = "🧬"
    formatted_message = f"{emoji} {message.upper()}"
    
    logging.info("")  # Optional line break before
    
    logging.info(f"{formatted_message}")
    logging.info("")

def print_stage_execution(message: str):
    emoji = "🚀"  
    formatted_message = f"{emoji}  {message.upper()}  {emoji}"
    separator = "█" * (len(formatted_message) + 6)

    logging.info("")  # Spacer before
    logging.info(separator)
    logging.info(f"█  {formatted_message}  █")
    logging.info(separator)
    logging.info("")  # Spacer after
