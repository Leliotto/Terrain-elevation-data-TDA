# Project Idea: TDA for Rescue Access Collapse in the Florida Keys

## Summary
Use the **Florida Keys** as a rescue-focused case study and ask:

**Can topological data analysis identify the inundation thresholds at which the Florida Keys’ dry land and evacuation road access fragment into isolated rescue zones, and can it locate the most critical access chokepoints?**

This is stronger than the military passability idea for your course because it is:
- more **structured and defensible**
- tied to **official public data**
- easier to validate visually and quantitatively
- a very natural fit for **persistence, lower-star filtrations, merge trees, and stability**

The project is not “just mapping.” It is a **topological analysis of accessibility collapse** as water rises.

## Key Changes
- **Case study and framing**
  - Use a **scenario-based** Florida Keys project, not a historical storm reconstruction.
  - Use Monroe County hurricane evacuation / rescue context as the operational motivation.
  - Treat the problem as: when does rising inundation break the Keys into disconnected dry-land / road-access components?

- **Main data**
  - **NOAA Sea Level Rise / coastal inundation data** for the Florida Keys, or equivalent NOAA Digital Coast elevation/inundation products.
  - **OpenStreetMap roads** for the evacuation spine and local access network.
  - Optional: **Monroe County evacuation information** as context, not as a primary analytic layer.
  - Optional: FEMA flood layers only if they help with discussion, not as a core dependency.

- **Main scalar/object to analyze**
  - Build a binary or weighted **accessibility surface** from:
    - dry land vs inundated land
    - road cells vs non-road cells
  - Recommended main object:
    - a **dry-road accessibility mask** at each inundation scenario
  - Optional richer version:
    - an **access cost surface** where inundated cells are forbidden, road cells are low-cost, and off-road dry terrain is higher-cost

- **TDA methods from the course**
  - **Lower-star / filtration viewpoint**: increase water level scenario by scenario and track how accessible space changes.
  - **`H0` persistence / connected components**: detect when one accessible region splits into multiple isolated rescue zones.
  - **Merge-tree style summary**: record when components disconnect or vanish as the scenario worsens.
  - **Stability analysis**: test whether the main chokepoints and breakup thresholds persist across resolution changes.
  - Optional extension if time permits:
    - `H1` on dry-land masks to study loops in accessible space
    - a Reeb-graph-like summary of the evacuation spine

- **Core outputs**
  - A map showing the **first major isolation thresholds**
  - A ranked list of **critical bridges / road segments / narrow connectors**
  - Persistence-style summaries of when communities or island clusters lose joint access
  - A stability result showing whether the main chokepoints remain under coarser resolution

## Implementation Changes
- **Study area**
  - Use a bounded AOI covering the main inhabited Florida Keys chain.
  - Keep one consistent AOI for all scenarios.

- **Scenario ladder**
  - Use fixed inundation levels, e.g. `0.5 ft`, `1 ft`, `2 ft`, `3 ft`, `4 ft`, `5 ft`, or the closest official NOAA increments available.
  - Each scenario produces:
    - dry-land mask
    - dry-road mask
    - connected-component summary
    - component sizes and isolation counts

- **Quantities to compute**
  - Number of accessible connected components (`β0`) per scenario
  - Size of largest accessible component
  - Number of communities / road clusters cut off from the main evacuation spine
  - First threshold where major fragmentation occurs
  - Persistence-like ranking of the most robust access corridors and the most critical failure points

- **Visuals to build**
  - Threshold montage of dry-road accessibility across scenarios
  - Component-count curve vs water level
  - Persistence / merge-tree summary for accessible regions
  - Geographic hotspot map of chokepoint failure locations
  - Zoom panels around the most critical bridges / narrow connectors

## Test Plan
- **Primary success criteria**
  - The method identifies a small number of visually plausible chokepoints rather than diffuse noise.
  - There is at least one clear fragmentation threshold where the main accessible region breaks apart.
  - The most important chokepoints remain broadly consistent under a coarser resolution rerun.

- **Robustness checks**
  - Rerun at a coarser grid or downsampled resolution.
  - Compare dry-land-only vs dry-land-plus-roads summaries.
  - Confirm that the main findings are not driven by one arbitrary threshold choice.

- **Interpretation checks**
  - The identified chokepoints should align with visible bridge/island connectors on the map.
  - The results should be explained as **accessibility fragmentation**, not as a full flood forecast.

## Assumptions and Defaults
- Use **Florida Keys** as the main case study.
- Use **scenario-based inundation**, not Hurricane Irma reconstruction.
- Use **terrain + roads** as the main object of analysis.
- Focus on **rescue access and evacuation fragmentation**, not building damage or hydrodynamic simulation.
- Keep `H0` / connectivity as the core analysis. Treat `H1`, Reeb-style summaries, and more advanced corridor skeletons as optional extensions if time remains.

## Public grounding
- Florida Keys evacuation / Monroe County emergency context: [Monroe County Emergency Management](https://www.monroecounty-fl.gov/907/Emergency-Management)
- NOAA sea-level / inundation planning context: [NOAA Sea Level Rise Viewer](https://coast.noaa.gov/slr/)
- NOAA Digital Coast data access: [Digital Coast](https://coast.noaa.gov/digitalcoast/)
- OpenStreetMap road data: [OpenStreetMap Export](https://www.openstreetmap.org/export)

## Recommended research question
**Can TDA detect the thresholds at which the Florida Keys’ dry-land and evacuation-road accessibility fragments into isolated rescue zones, and can it identify the most critical chokepoints in that process?**
