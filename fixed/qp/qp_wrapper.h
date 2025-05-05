#pragma once

#include <tuple>
#include <Eigen/Sparse>
#include <Eigen/Dense>


using Eigen::SparseMatrix;

/// @brief Преобразует разреженную матрицу Eigen в формат CCS
/// Формат CCS (он же CSC) соответствует ColMajor
/// Если матрица RowMajor, надо ее сначала преобразовать в ColMajor
/// https://eigen.tuxfamily.org/dox/group__TutorialSparse.html (раздел  Sparse matrix format)
/// @param Матрица
/// @return [values; rows; col_indices]
/// values - ненулевые элементы
/// rows - индексы их строк
/// col_indices - индексы начала столбцов в values и rows (отличается по размерности от этих массивов)
/// Пример из Deepseek
/// [2  0  3]
/// [0  5  0]
/// [7  0  9]
/// values[] = {2, 7, 5, 3, 9} 
/// rows[]   = {0, 2, 1, 0, 2} 
/// col_indices[] = {0, 2, 3, 5} 
template <typename IndexType = long long /*aka qp_int*/>
inline std::tuple<
    std::vector<double>, // aka qp_real
    std::vector<IndexType>,
    std::vector<IndexType>>
    get_sparse_matrix_CCS(const SparseMatrix<double, Eigen::ColMajor>& _matrix)
{
    SparseMatrix<double, Eigen::ColMajor> matrix = _matrix;
    matrix.makeCompressed();

    std::vector<double> values;
    std::vector<IndexType> rows;
    std::vector<IndexType> cols;

    for (IndexType k = 0; k < matrix.outerSize(); ++k) {
        typename SparseMatrix<double, Eigen::ColMajor>::InnerIterator it(matrix, k);

        cols.push_back(values.size());

        for (; it; ++it) {
            values.push_back(it.value());
            rows.push_back(static_cast<IndexType>(it.row()));
        }
    }
    cols.push_back(values.size());

    return std::make_tuple(std::move(values), std::move(rows), std::move(cols));
}


#ifdef FIXED_CONSTRAINED_EQUATION_SOLVER
template <typename IndexType>
VectorXd solve_quadprog_inequalities_swift(
    const SparseMatrix<double, ColMajor, IndexType>& H, const VectorXd& f,
    const SparseMatrix<double, ColMajor, IndexType>& A, const VectorXd& b);
#endif
