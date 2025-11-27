# Single Cell Browser Rework

Aim is to rework the existing single cell browser in a way where it can be easily extended and therefore not require 
costly changes for each experiment.

## Architecture

The architecture of the app will follow a hexagonal pattern with the data objects and business logic within the `core`.
Outside of that will be the service layer and then the adapters and ports to real world applications.

### Key Abstractions

The primary use case for abstraction in this app is for the figures, `view_base.py` provides a set of methods and properties
for subtypes to follow. This enforces behaviour across subtypes as well as making it easier to implement new ones in the future.

Also, `view_registry.py` 