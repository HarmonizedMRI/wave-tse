# Wave TSE

Pulseq-based 2D Wave-CAIPI turbo spin echo (TSE) sequence generation for
single-band (MB1) and two-band simultaneous-multislice (MB2) acquisitions,
a mixed-dimensional FLASH calibration, and Python reconstruction from
Siemens TWIX data.

The repository has two intentionally separate parts:

- MATLAB generates the MB2 TSE, the MB1 TSE and the calibration as three
  independent Pulseq sequences.
- Python reconstructs TSE data from the calibration TWIX, the TSE TWIX, and
  the matching `.seq` files and sidecars.

The repository is being assembled; the sequence and reconstruction guides
will be added under `docs/`.

## License

MIT License. See `LICENSE`.
