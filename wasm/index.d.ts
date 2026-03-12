/**
 * HPOP SDN Plugin — TypeScript declarations
 *
 * High-Precision Orbit Propagator compiled to WebAssembly.
 * Parses VCM text/binary, propagates orbits, outputs OEM/OCM FlatBuffers.
 */

interface HPOPModule {
    /**
     * Propagate orbit from JSON input → OEM FlatBuffer binary string.
     * JSON fields: epoch, x,y,z, vx,vy,vz, duration_days, step_seconds,
     *              mass_kg, drag_area, cd, srp_area, cr, f107
     */
    propagate(json: string): string;

    /**
     * Propagate from VCM-schema JSON → OEM FlatBuffer binary string.
     * Accepts both simplified and VCM-schema field names.
     */
    propagateVCM(json: string): string;

    /**
     * Parse raw 18 SPCS text VCM → OCM FlatBuffer binary (Uint8Array).
     * Input: VCM text with <> line prefixes
     * Output: $OCM FlatBuffer binary
     */
    parseVCM(vcmText: string): Uint8Array;

    /**
     * Convert OCM FlatBuffer binary → raw text VCM string.
     * Input: Uint8Array of $OCM FlatBuffer
     * Output: 18 SPCS ASCII VCM text
     */
    convertOCM(ocmBinary: Uint8Array): string;

    /** Get plugin version string */
    getVersion(): string;

    /** Get available force model names as JSON array string */
    getForceModels(): string;

    /** Allocate memory on WASM heap */
    _malloc(size: number): number;

    /** Free WASM heap memory */
    _free(ptr: number): void;

    /** WASM heap views */
    HEAPU8: Uint8Array;
    HEAPF64: Float64Array;
}

/**
 * Initialize the HPOP WASM module.
 * @returns Promise that resolves to the initialized module.
 */
declare function HPOPWasmInit(): Promise<HPOPModule>;
export default HPOPWasmInit;
export { HPOPModule };
