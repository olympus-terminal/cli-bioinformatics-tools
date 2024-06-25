#!/bin/bash

# Array of Perl modules to install
modules=(
    "Text::Format"
    "Statistics::Basic"
    "Readonly"
    "Set::IntSpan"
    "Font::TTF"
    "Config::General"
    "GD::Polyline"
    "GD"
    "Math::VecStat"
    "Math::Bezier"
    "SVG"
)

# Function to install a module
install_module() {
    echo "Installing $1..."
    cpanm "$1"
    if [ $? -eq 0 ]; then
        echo "$1 installed successfully."
    else
        echo "Failed to install $1. Please check the error and try manually."
    fi
    echo
}

# Main installation loop
for module in "${modules[@]}"; do
    install_module "$module"
done

echo "Installation process completed."
