#ifndef AI_H
#define AI_H

#include <iostream>
#include <limits>
#include <vector>
#include "Board.h"
#include "Ply.h"
#include "Rectangle.h"

class AI {
public:
	AI() = delete;
	static bool gameWon(Board *board, std::vector<int> &bitBoard);
	static bool gameWon(Board *board, bool isBlack);
	static Position computeScore(Board *board, int8_t maxSearchDepth, double maxSearchTime, bool isBlack, bool printAIInformation);
	static int16_t evaluateBoard(Board *board, bool isBlack, bool printInformation);
private:
	static void computeNextPositions(Board *board, std::vector<std::vector<Position>> &nextPositions, int hashIndexPrimary, int hashIndexSecondary, int8_t maxLevels, int8_t levelsToBottom, bool isBlack);
	static int setupHashObject(Board *board, int hashIndexPrimary);
	static int16_t computeScoreUsingAlphaBeta(Board *board, std::vector<std::vector<Position>> &nextPositions, int8_t maxLevels, int8_t levelsToBottom, int16_t alpha, int16_t beta, bool black, bool mostComputeScore);
	static void alphaBetaInnerLoop(Board *board, std::vector<std::vector<Position>> &nextPositions, Position position, int16_t *alpha, int16_t *beta, int8_t levelsToBottom, int8_t maxLevels, bool isBlack, bool mostComputeScore, Ply *bestNextPly, bool *hasFoundABestNextPly);
	static int numberOfFoursOpen(Board *board, std::vector<int> &bitBoard, std::vector<int> &bitBoardEmpty);
	static bool setwinningPosition(int leftColumnLimit, int rightColumnLimit, int row, int columnNumber, int columnModifier, Position &winningPosition);
	static bool hasAWinningPly(Board *board, std::vector<int> &bitBoard, std::vector<int> &bitBoardEmpty, Position &winningPosition);
	static int numberOfFoursHalfOpen(Board *board, std::vector<int> &bitBoard, std::vector<int> &bitBoardEmpty);
	static int numberOfThreesOpen(Board *board, std::vector<int> &bitBoard, std::vector<int> &bitBoardEmpty);
	static int numberOfThreesHalfOpen(Board *board, std::vector<int> &bitBoard, std::vector<int> &bitBoardEmpty);
	static int numberOfTwosOpen(Board *board, std::vector<int> &bitBoard, std::vector<int> &bitBoardEmpty);
	static int numberOfTwosHalfOpen(Board *board, std::vector<int> &bitBoard, std::vector<int> &bitBoardEmpty);
};

#endif
