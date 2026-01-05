# Ketcher Function Guide for AI Tool Calling

This guide provides an overview of functions within the Ketcher codebase (specifically `ketcher-core`) that are suitable for exposure as tools for an AI agent. These functions allow for programmatic manipulation and analysis of chemical structures.

## Overview

The core logic of Ketcher resides in the `ketcher-core` package. The primary entry point for interaction is the `Ketcher` class, which coordinates the editor, struct service (Indigo), and other components.

Most chemistry-specific operations (aromatization, layout, calculation) are handled by the `Indigo` service wrapper, which communicates with the backend (or WASM) service.

## Key Classes

*   **`Ketcher`**: Located in `ketcher-core/src/application/ketcher.ts`. This is the main controller. It provides high-level methods to get/set molecules, export images, and access the `Indigo` service.
*   **`Indigo`**: Located in `ketcher-core/src/application/indigo.ts`. This class wraps the `StructService` and provides specific chemical operations.

## Potential AI Tools

The following functions are good candidates for AI tools. They generally accept a molecule structure (often in KET or SMILES format) and return a modified structure or calculated property.

### 1. Aromatize Structure

Converts a molecule to its aromatic form (e.g., benzene ring with circle or aromatic bonds).

*   **Function**: `ketcher.indigo.aromatize(struct)`
*   **Description**: Aromatizes the given chemical structure.
*   **Parameters**:
    *   `struct` (string | Struct): The chemical structure. Can be a string (SMILES, Molfile, etc.) or a Ketcher `Struct` object.
*   **Return Type**: `Promise<Struct>` (The aromatized structure)
*   **Source**: `packages/ketcher-core/src/application/indigo.ts`

### 2. Dearomatize Structure

Converts a molecule to its Kekulé form (alternating single and double bonds).

*   **Function**: `ketcher.indigo.dearomatize(struct)`
*   **Description**: Dearomatizes the given chemical structure.
*   **Parameters**:
    *   `struct` (string | Struct): The chemical structure.
*   **Return Type**: `Promise<Struct>`
*   **Source**: `packages/ketcher-core/src/application/indigo.ts`

### 3. Auto Layout (Clean Up)

Automatically arranges the atoms and bonds of a molecule in a readable 2D layout.

*   **Function**: `ketcher.indigo.layout(struct)` or `ketcher.layout()` (which acts on current editor content)
*   **Description**: Performs 2D layout on the structure.
*   **Parameters**:
    *   `struct` (string | Struct): The chemical structure.
*   **Return Type**: `Promise<Struct>`
*   **Source**: `packages/ketcher-core/src/application/indigo.ts`

### 4. Calculate Properties

Calculates various chemical properties (Molecular Weight, Mass, Formula, etc.).

*   **Function**: `ketcher.indigo.calculate(struct, options)`
*   **Description**: Calculates specified properties for the structure.
*   **Parameters**:
    *   `struct` (string | Struct): The chemical structure.
    *   `options` (object): Optional. `{ properties: string[] }`. Default properties include 'molecular-weight', 'most-abundant-mass', 'monoisotopic-mass', 'gross', 'mass-composition'.
*   **Return Type**: `Promise<CalculateResult>` (Object with calculated values)
*   **Source**: `packages/ketcher-core/src/application/indigo.ts`

### 5. Check Structure

Checks the structure for errors or warnings (valence, radicals, etc.).

*   **Function**: `ketcher.indigo.check(struct, options)`
*   **Description**: Validates the structure against a set of rules.
*   **Parameters**:
    *   `struct` (string | Struct): The chemical structure.
    *   `options` (object): Optional. `{ types: string[] }`. Types can be 'valence', 'radicals', 'pseudoatoms', etc.
*   **Return Type**: `Promise<CheckResult>` (Map of errors/warnings)
*   **Source**: `packages/ketcher-core/src/application/indigo.ts`

### 6. Get Structure (Format Conversion)

Retrieves the current structure in a specific format (SMILES, Molfile, InChI, etc.).

*   **Function**: `ketcher.getSmiles()`, `ketcher.getMolfile()`, `ketcher.getInchi()`, etc.
*   **Description**: Exports the current editor structure to a string in the requested format.
*   **Parameters**:
    *   Various, depending on the format (e.g., `isExtended` for SMILES).
*   **Return Type**: `Promise<string>`
*   **Source**: `packages/ketcher-core/src/application/ketcher.ts`

### 7. Set Molecule

Loads a molecule into the Ketcher editor.

*   **Function**: `ketcher.setMolecule(structStr)`
*   **Description**: Parses and renders the given structure string in the editor.
*   **Parameters**:
    *   `structStr` (string): The structure string (SMILES, Molfile, etc.).
*   **Return Type**: `Promise<void>`
*   **Source**: `packages/ketcher-core/src/application/ketcher.ts`

### 8. Automap Reaction

Automatically maps atoms in a chemical reaction.

*   **Function**: `ketcher.indigo.automap(struct, options)`
*   **Description**: automatically maps atoms in a chemical reaction.
*   **Parameters**:
    *   `struct` (string | Struct): The reaction structure.
    *   `options` (object): `{ mode: 'discard' | 'keep' | 'alter' | 'clear' }`
*   **Return Type**: `Promise<Struct>`
*   **Source**: `packages/ketcher-core/src/application/indigo.ts`

### 9. Calculate CIP (Stereochemistry)

Calculates Cahn-Ingold-Prelog stereochemistry descriptors (R/S, E/Z).

*   **Function**: `ketcher.indigo.calculateCip(struct)`
*   **Description**: Calculates CIP descriptors.
*   **Parameters**:
    *   `struct` (string | Struct): The chemical structure.
*   **Return Type**: `Promise<Struct>`
*   **Source**: `packages/ketcher-core/src/application/indigo.ts`

## Implementation Guide for AI Agent

To implement these as tools for an AI agent:

1.  **Access**: Ensure the agent has access to the `ketcher` instance (usually available globally as `window.ketcher` in the browser environment if exposed, or via the React context).
2.  **Wrappers**: Create wrapper functions that bridge the AI's string-based tool calling with Ketcher's object-based API.
    *   *Example Wrapper*:
        ```javascript
        async function aromatizeTool(smiles) {
            // 1. Load into Ketcher (optional, if you want to visualize)
            // await ketcher.setMolecule(smiles);

            // 2. Perform operation
            // Note: ketcher.indigo.aromatize takes a structure.
            // If passing a string, it might need to be in a format Ketcher understands or the internal Struct object.
            // The indigo.ts wrapper handles string inputs by converting them if needed.
            const resultStruct = await ketcher.indigo.aromatize(smiles);

            // 3. Convert back to string (e.g., SMILES) to return to AI
            // We might need a helper to convert Struct -> SMILES if not directly supported by the return of aromatize
            // ketcher.indigo.aromatize returns a Struct object (which can be serialized).
            // To get SMILES, you might need to load it into ketcher and call getSmiles(),
            // or use a separate conversion service if available.

            // Simpler flow if visualizing:
            await ketcher.setMolecule(smiles);
            const aromatized = await ketcher.indigo.aromatize(ketcher.editor.struct());
            ketcher.editor.struct(aromatized); // Update editor
            return await ketcher.getSmiles();
        }
        ```
3.  **Context**: The agent should know that these operations might modify the state of the editor if used directly on `ketcher.editor.struct()`.

## File Index

*   `packages/ketcher-core/src/application/ketcher.ts`: Main API surface.
*   `packages/ketcher-core/src/application/indigo.ts`: Chemistry logic wrapper.
*   `packages/ketcher-core/src/domain/services/struct/structService.types.ts`: Type definitions for inputs/outputs.
