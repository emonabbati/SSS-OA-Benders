# Resilient Supplier Selection and Order Allocation

Implementation of a three-stage framework for resilient supplier selection under disruption risks, using:

- Bayesian Networks for scenario generation
- Markov Chains for availability modeling
- Bi-objective optimization with Benders Decomposition

## Requirements

- Julia 1.11+
- CPLEX 12.10+ (with JuMP interface)

## Results

The code generates:
- Pareto frontier (Figure 5)
- Availability curves (Figure 6)
- Robust vs. deterministic comparison (Figures 7-8)
- Sensitivity analysis (Figures 9-12)

## Citation

If you use this code in your research, please cite:

```bibtex
@article{alipour2026resilient,
  title={Resilient supplier selection and order allocation problem integrated by Markov chain and Bayesian network-based disruption modeling},
  author={Alipour, Z. and Monabbati, S.E.},
  note={Under review},
  year={2026}
}
```

## License

MIT License - see the [LICENSE](LICENSE) file for details.